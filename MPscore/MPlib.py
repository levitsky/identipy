import sys
import re
import SSRCalc
from string import punctuation, lowercase
from pyteomics import parser, biolccc, pepxml
from pyteomics import mass
from pyteomics.auxiliary import linear_regression
from pyteomics import achrom
import numpy as np
from scipy.stats import scoreatpercentile
from copy import copy
from scipy.spatial import cKDTree
from collections import Counter

try:
    from configparser import RawConfigParser
except ImportError:
    from ConfigParser import RawConfigParser

def get_dbname(prot, pepxml_type='tandem'):
    if pepxml_type != 'omssa':
        try:
            if not any(prot['protein'].startswith(tag) for tag in ['sp', 'tr', 'DECOY_sp', 'DECOY_tr']):
                if any(prot['protein_descr'].startswith(tag) for tag in ['SWISS-PROT:', 'TREMBL:']):
                    return prot['protein']
                if '|' not in prot['protein']:
                    return str(prot['protein']+' '+prot['protein_descr']).replace('DECOY_', '')
                else:
                    return str(prot['protein']+'>'+prot['protein_descr']).replace('DECOY_', '')
            else:
                return str(prot['protein'].split('|')[1]).replace('DECOY_', '')
        except:
            return str(prot['protein'].split('_')[0]).replace('DECOY_', '')
    else:
        return str(prot['protein_descr'].split('|')[1]).replace('DECOY_', '')


def get_aa_mass(settings):
    aa_mass = mass.std_aa_mass.copy()
    fmods = settings.get('modifications', 'fixed')
    if fmods:
        for mod in re.split(r'[,;]\s*', fmods):
            m, aa = parser._split_label(mod)
            aa_mass[aa] += settings.getfloat('modifications', m)
            aa_mass[mod] = aa_mass[aa] + settings.getfloat('modifications', m)
    vmods = settings.get('modifications', 'variable')
    if vmods:
        mods = [(l[:-1], l[-1]) for l in re.split(r',\s*', vmods)]
        for (mod, aa), char in zip(mods, punctuation):
            if aa == '[':
                aa_mass[mod + '-'] = settings.getfloat('modifications', mod) + settings.getfloat('modifications', 'protein nterm cleavage')
            elif aa == ']':
                aa_mass['-' + mod] = settings.getfloat('modifications', mod) + settings.getfloat('modifications', 'protein cterm cleavage')
            else:
                aa_mass[mod + aa] = aa_mass[aa] + settings.getfloat('modifications', mod)
    return aa_mass

def filter_evalue_prots(prots, FDR=1.0):
    target_evalues = np.array([v['expect'] for k, v in prots.iteritems() if not k.startswith('L')])
    decoy_evalues = np.array([v['expect'] for k, v in prots.iteritems() if k.startswith('L')])
    target_evalues.sort()
    best_cut_evalue = None
    real_FDR = 0
    for cut_evalue in target_evalues:
        counter_target = target_evalues[target_evalues <= cut_evalue].size
        counter_decoy = decoy_evalues[decoy_evalues <= cut_evalue].size
        if counter_target and (float(counter_decoy) / float(counter_target)) * 100 <= float(FDR):
            best_cut_evalue = cut_evalue
            real_FDR = round(float(counter_decoy) / float(counter_target) * 100, 1)
    if not best_cut_evalue:
        best_cut_evalue = 0
    print real_FDR, best_cut_evalue, 'protein e-value'
    new_prots = {}
    for k, v in prots.iteritems():
        if v['expect'] <= best_cut_evalue and not k.startswith('L'):
            new_prots[k] = v
    return new_prots

def get_settings(fname=None, default_name='default.cfg'):
    """Read a configuration file and return a :py:class:`RawConfigParser` object.
    """
    kw = {'inline_comment_prefixes': ('#', ';')
            } if sys.version_info.major == 3 else {}
    kw['dict_type'] = dict
    kw['allow_no_value'] = True

    raw_config = RawConfigParser(**kw)
    if default_name:
        raw_config.read(default_name)
    if fname:
        raw_config.read(fname)
    return raw_config


def FDbinSize(X):
    """Calculates the Freedman-Diaconis bin size for
    a data set for use in making a histogram
    Arguments:
    X:  1D Data set
    Returns:
    h:  F-D bin size
    """
    X = np.sort(X)
    upperQuartile = scoreatpercentile(X, 75)
    lowerQuartile = scoreatpercentile(X, 25)
    IQR = upperQuartile - lowerQuartile
    h = 2. * IQR / len(X) ** (1. / 3.)
    return h

class Descriptor():
    def __init__(self, name, formula, group='A', binsize='auto'):
        self.name = name
        self.formula = formula
        self.group = group
        self.binsize = binsize
        self.normK = None
        self.array = None

    def get_array(self, peptides):
        tmp = np.array([self.formula(peptide) for peptide in peptides.peptideslist])
        tmp.sort()
        self.array = tmp
        return tmp

    def get_binsize(self, peptides=None):
        if self.binsize == 'auto':
            if self.group == 'B':
                return 0.25
            else:
                return FDbinSize(self.get_array(peptides))
        else:
            return self.binsize


class PeptideList:
    def __init__(self, settings=None, mods=None):
        self.peptideslist = []
        self.calibrate_coeff = None
        self.pepxml_type = ''
        self.RC = False
        self.settings = settings
        if not mods:
            self.modification_list = {}
        else:
            self.modification_list = mods
        self.total_number_of_PSMs = 0
        self.total_number_of_PSMs_decoy = 0
        self.total_number_of_peptides_in_searchspace = 0
        self.total_number_of_proteins_in_searchspace = 0
        self.total_number_of_spectra = 0
        self.nterm_mass = self.settings.getfloat('modifications', 'protein nterm cleavage')
        self.cterm_mass = self.settings.getfloat('modifications', 'protein cterm cleavage')
        self.aa_list = get_aa_mass(settings)
        fmods = self.settings.get('modifications', 'fixed')
        if fmods:
            for mod in re.split(r'[,;]\s*', fmods):
                m, aa = parser._split_label(mod)
                self.modification_list[str(int(mass.std_aa_mass[aa] + settings.getfloat('modifications', m)))] = m
        vmods = settings.get('modifications', 'variable')
        if vmods:
            mods = [(l[:-1], l[-1]) for l in re.split(r',\s*', vmods)]
            for (mod, aa), char in zip(mods, punctuation):
                if aa == '[':
                    self.modification_list[str(int(self.nterm_mass + settings.getfloat('modifications', mod)))] = mod + '-'
                elif aa == ']':
                    self.modification_list[str(int(self.cterm_mass + settings.getfloat('modifications', mod)))] = '-' + mod
                else:
                    self.modification_list[str(int(mass.std_aa_mass[aa] + settings.getfloat('modifications', mod)))] = mod

    def __len__(self):
        return len(self.peptideslist)

    def get_number_of_peptides(self):
        return len(set(p.sequence for p in self.peptideslist))

    def get_number_of_spectra(self):
        """Returns the number of MS/MS spectra used for the search. If mgf file is not available,
         returns number of identified PSMs as approximation.
        """
        if self.total_number_of_spectra:
            return self.total_number_of_spectra
        else:
            return self.total_number_of_PSMs

    def get_from_pepxmlfile(self, pepxmlfile, min_charge=1, max_charge=0):
        for line in open(pepxmlfile, 'r'):
            if line.startswith('<search_summary') or line.startswith('    <search_summary'):
                if "X! Tandem" in line:
                    self.pepxml_type = 'tandem'
                elif "OMSSA" in line:
                    self.pepxml_type = 'omssa'
                elif "MASCOT" in line:
                    self.pepxml_type = 'mascot'
                else:
                    print "Unknown search_engine"
                break

        pepxml_params = {k: v for d in pepxml.iterfind(pepxmlfile, 'parameter name') for k, v in d.items()}
        self.total_number_of_peptides_in_searchspace = int(pepxml_params.get('modelling, total peptides used', self.total_number_of_peptides_in_searchspace))
        self.total_number_of_proteins_in_searchspace = int(pepxml_params.get('modelling, total proteins used', self.total_number_of_proteins_in_searchspace))
        self.total_number_of_spectra = int(pepxml_params.get('modelling, total spectra used', self.total_number_of_spectra))

        standard_aminoacids = set(k for k in mass.std_aa_comp if '-' not in k)
        first_psm = True
        for record in pepxml.read(pepxmlfile):
            if 'search_hit' in record:
                if int(min_charge) <= int(record['assumed_charge']) and (int(record['assumed_charge']) <= int(max_charge) or not max_charge):
                    if first_psm:
                        if 'num_missed_cleavages' not in record['search_hit'][0]:
                            print 'missed cleavages are missed in pepxml file, using 0 value for all peptides'
                        try:
                            float(record['retention_time_sec'])
                        except:
                            try:
                                float(record['spectrum'].split(',')[2].split()[0])
                            except:
                                print 'RT experimental is missed in pepxml file, using 0 value for all peptides'
                        first_psm = False

                    sequence = record['search_hit'][0]['peptide']
                    if not set(sequence).difference(standard_aminoacids):
                        mc = record['search_hit'][0].get('num_missed_cleavages', 0)
                        modified_code = record['search_hit'][0]['modified_peptide']
                        modifications = record['search_hit'][0]['modifications']
                        try:
                            evalue = record['search_hit'][0]['search_score']['expect']
                        except:
                            try:
                                evalue = 1.0 / float(record['search_hit'][0]['search_score']['ionscore'])
                            except IOError:
                                'Cannot read e-value!'
                        try:
                            sumI = round(10 ** float(record['search_hit'][0]['search_score']['sumI']), 0)
                        except:
                            sumI = 0
                        spectrum = record['spectrum']
                        pcharge = record['assumed_charge']
                        mass_exp = record['precursor_neutral_mass']

                        pept = Peptide(sequence=sequence, settings=self.settings, modified_code=modified_code, evalue=evalue, spectrum=spectrum, pcharge=pcharge, mass_exp=mass_exp, modifications=modifications, modification_list=self.modification_list, custom_aa_mass=self.aa_list, sumI=sumI, mc=mc)
                        try:
                            pept.RT_exp = float(record['retention_time_sec']) / 60
                        except:
                            try:
                                pept.RT_exp = float(record['spectrum'].split(',')[2].split()[0])
                            except:
                                pept.RT_exp = 0

                        decoy_tags = [':reversed', 'DECOY_', 'rev_', 'Random sequence.']
                        if any([all([all((not protein.get(key, '') or not protein[key].startswith(tag)) and (not protein.get(key, '') or not protein[key].endswith(tag)) for key in ['protein', 'protein_descr']) for tag in decoy_tags]) for protein in record['search_hit'][0]['proteins']]):
                            pept.note = 'target'
                        else:
                            pept.note = 'decoy'

                        for prot in record['search_hit'][0]['proteins']:
                            if get_dbname(prot, self.pepxml_type) not in [protein.dbname for protein in pept.parentproteins]:
                                pept.parentproteins.append(Protein(dbname=get_dbname(prot, self.pepxml_type), description=prot.get('protein_descr', None)))

                        if len(pept.parentproteins) and (not modifications or Counter(v['position'] for v in modifications).most_common(1)[0][1] <= 2):
                            self.peptideslist.append(pept)


    def modified_peptides(self):
        for peptide in self.peptideslist:
            peptide.modified_peptide()

    def get_RC(self):
        seqs = [pept.modified_sequence for pept in self.peptideslist]
        RTexp = [pept.RT_exp for pept in self.peptideslist]
        RC_def = achrom.RCs_gilar_rp
        aa_labels = set(RC_def['aa'].keys())
        for pept in self.peptideslist:
            for v in pept.modification_list.itervalues():
                aa_labels.add(v)
        xdict = {}
        for key, val in RC_def['aa'].items():
            xdict[key] = [val, None]
        RC_dict = achrom.get_RCs_vary_lcp(seqs, RTexp, labels=aa_labels)
        for key, val in RC_dict['aa'].items():
            try:
                xdict[key][1] = val
            except:
                xdict[key] = [None, val]
        a, b, _, _ = linear_regression([x[0] for x in xdict.values() if all(v != None for v in x)], [x[1] for x in xdict.values() if all(v != None for v in x)])
        for key, x in xdict.items():
            if x[1] == None:
                x[1] = x[0] * a + b
            RC_dict['aa'][key] = x[1]
        if 'C' not in RC_dict['aa']:
            RC_dict['aa']['C'] = RC_dict['aa']['C*']
        self.RC = RC_dict

    def calc_RT(self, calibrate_coeff=(1, 0, 0, 0), RTtype='achrom'):
        if RTtype == 'ssrcalc':
            ps = list(set([peptide.sequence for peptide in self.peptideslist]))
            SSRCalc_RTs = SSRCalc.calculate_RH(ps[:], pore_size=100, ion_pairing_agent='FA')

        for peptide in self.peptideslist:
            if RTtype == 'achrom':
                peptide.RT_predicted = achrom.calculate_RT(peptide.modified_sequence, self.RC, raise_no_mod=False)
                if np.isinf(peptide.RT_predicted):
                    elems = peptide.modified_sequence.split('-')
                    if len(elems) > 1:
                        if not all(el in parser.std_amino_acids for el in elems[0]) and str(elems[0] + '-') not in self.RC['aa']:
                            self.RC['aa'][str(elems[0] + '-')] = self.RC['aa']['H-']
                        elif not all(el in parser.std_amino_acids for el in elems[-1]) and str('-' + elems[-1]) not in self.RC['aa']:
                            self.RC['aa'][str('-' + elems[-1])] = self.RC['aa']['-OH']
                    peptide.RT_predicted = achrom.calculate_RT(peptide.modified_sequence, self.RC, raise_no_mod=False)

            elif RTtype == 'ssrcalc':
                SSRCalc_RT = SSRCalc_RTs[peptide.sequence]
                if SSRCalc_RT is not None:
                    peptide.RT_predicted = float(SSRCalc_RT) * calibrate_coeff[0] + calibrate_coeff[1]
                else:
                    peptide.RT_predicted = 0
                    print 'SSRCalc error'
            elif RTtype == 'biolccc':
                peptide.RT_predicted = biolccc.calculateRT(peptide.sequence, biolccc.rpAcnTfaChain, biolccc.standardChromoConditions)
            else:
                print 'RT_type error'

    def filter_RT(self, RT_tolerance):
        j = len(self.peptideslist) - 1
        while j >= 0:
            if abs(float(self.peptideslist[j].RT_predicted) - float(self.peptideslist[j].RT_exp)) > float(RT_tolerance):
                self.peptideslist.pop(j)
            j -= 1

    def get_calibrate_coeff(self):
        peptides = []
        peptides_added = {}
        for peptide in self.peptideslist:
            if peptide.sequence not in peptides_added:
                peptides_added[peptide.sequence] = [peptide.RT_exp, ]
                peptides.append([peptide.RT_predicted, peptide.RT_exp])
            else:
                if any(abs(peptide.RT_exp - v) < 2 for v in peptides_added[peptide.sequence]):
                    pass
                else:
                    peptides_added[peptide.sequence].append(peptide.RT_exp)
                    peptides.append([peptide.RT_predicted, peptide.RT_exp])
        aux_RT = linear_regression([val[0] for val in peptides], [val[1] for val in peptides])
        return aux_RT

    def filter_modifications(self, RT_type=None):
        j = len(self.peptideslist) - 1
        while j >= 0:
            if self.peptideslist[j].modified_code.count('[') - sum(self.peptideslist[j].modified_code.count('[%s]' % (x, )) for x in ([160, 181, 166, 243] if RT_type=='biolccc' else [160, ]) ) != 0:
                self.peptideslist.pop(j)
            j -= 1

    def filter_decoy(self):
        j = len(self.peptideslist) - 1
        while j >= 0:
            if self.peptideslist[j].note == 'decoy':
                self.peptideslist.pop(j)
            j -= 1

    def filter_evalue_new(self, FDR=1, FDR2=1, useMP=True, drop_decoy=True, toprint=False):
        "A function for filtering PSMs by e-value and MP-score with some FDR"
        target_evalues, decoy_evalues = [], []
        for peptide in self.peptideslist:
            if peptide.note == 'target':
                target_evalues.append(float(peptide.evalue))
            elif peptide.note == 'decoy':
                decoy_evalues.append(float(peptide.evalue))
        target_evalues = np.array(target_evalues)
        decoy_evalues = np.array(decoy_evalues)
        target_evalues.sort()
        best_cut_evalue = None
        real_FDR = 0
        for cut_evalue in target_evalues:
            counter_target = target_evalues[target_evalues <= cut_evalue].size
            counter_decoy = decoy_evalues[decoy_evalues <= cut_evalue].size
            if counter_target and (float(counter_decoy) / float(counter_target)) * 100 <= float(FDR):
                best_cut_evalue = cut_evalue
                real_FDR = round(float(counter_decoy) / float(counter_target) * 100, 1)
        if not best_cut_evalue:
            best_cut_evalue = 0
        print real_FDR, best_cut_evalue, 'e-value'

        best_cut_peptscore = 1.1
        if useMP:
            target_peptscores, decoy_peptscores = [], []
            for peptide in self.peptideslist:
                if peptide.evalue >= best_cut_evalue:
                    if peptide.note == 'target':
                        target_peptscores.append(float(peptide.peptscore))
                    elif peptide.note == 'decoy':
                        decoy_peptscores.append(float(peptide.peptscore))
            target_peptscores = np.array(target_peptscores)
            decoy_peptscores = np.array(decoy_peptscores)
            target_peptscores = np.sort(target_peptscores)[::-1]
            real_FDR = 0
            for cut_peptscore in target_peptscores:
                counter_target = target_peptscores[target_peptscores >= cut_peptscore].size
                counter_decoy = decoy_peptscores[decoy_peptscores >= cut_peptscore].size
                if counter_target and (float(counter_decoy) / float(counter_target)) * 100 <= float(FDR2):
                    best_cut_peptscore = cut_peptscore
                    real_FDR = round(float(counter_decoy) / float(counter_target) * 100, 1)
            print real_FDR, best_cut_peptscore, 'MP score'
        new_peptides = self.copy_empty()
        for peptide in self.peptideslist:
            if peptide.evalue <= best_cut_evalue or (useMP and peptide.peptscore >= best_cut_peptscore):
                new_peptides.peptideslist.append(peptide)
        if drop_decoy:
            new_peptides.filter_decoy()
        return (new_peptides, best_cut_evalue, best_cut_peptscore)

    def copy_empty(self):
        new_peptides = PeptideList(self.settings)
        new_peptides.pepxml_type = self.pepxml_type
        new_peptides.total_number_of_spectra = self.total_number_of_spectra
        new_peptides.total_number_of_PSMs = self.total_number_of_PSMs
        new_peptides.total_number_of_PSMs_decoy = self.total_number_of_PSMs_decoy
        new_peptides.total_number_of_proteins_in_searchspace = self.total_number_of_proteins_in_searchspace
        new_peptides.total_number_of_peptides_in_searchspace = self.total_number_of_peptides_in_searchspace
        new_peptides.calibrate_coeff = self.calibrate_coeff
        new_peptides.RC = self.RC
        new_peptides.modification_list = self.modification_list
        return new_peptides

    def update(self, new_peptides):
        self.pepxml_type = new_peptides.pepxml_type
        self.settings = new_peptides.settings
        self.pepxml_type = new_peptides.pepxml_type
        self.total_number_of_spectra = new_peptides.total_number_of_spectra
        self.total_number_of_PSMs = new_peptides.total_number_of_PSMs
        self.total_number_of_PSMs_decoy = new_peptides.total_number_of_PSMs_decoy
        self.total_number_of_proteins_in_searchspace = max(new_peptides.total_number_of_proteins_in_searchspace, self.total_number_of_proteins_in_searchspace)
        self.total_number_of_peptides_in_searchspace = max(new_peptides.total_number_of_peptides_in_searchspace, self.total_number_of_peptides_in_searchspace)
        self.calibrate_coeff = new_peptides.calibrate_coeff
        self.RC = new_peptides.RC
        self.modification_list = new_peptides.modification_list
        self.peptideslist.extend(new_peptides.peptideslist)

    def remove_duplicate_sequences(self):
        new_peptides = self.copy_empty()
        for peptide in self.peptideslist:
            if peptide.sequence not in set(p.sequence for p in new_peptides.peptideslist) and peptide.evalue == min([p.evalue for p in self.peptideslist if p.sequence == peptide.sequence]):
                new_peptides.peptideslist.append(peptide)
        return new_peptides

class Protein:
    def __init__(self, dbname, description='Unknown'):
        self.dbname = dbname
        self.description = description

class Peptide:
    def __init__(self, sequence, settings, modified_code='', pcharge=0, RT_exp=False, evalue=0, note='unknown', spectrum='', mass_exp=0, modifications=[], modification_list={}, custom_aa_mass=None, sumI=0, mc=None):
        self.sequence = sequence
        self.modified_code = modified_code
        self.modified_sequence = sequence
        self.modifications = modifications
        self.modification_list = modification_list
        self.pcharge = int(pcharge)
        self.nomodifications = 0
        self.aa_mass = custom_aa_mass
        self.pmass = float(mass.calculate_mass(sequence=self.sequence, charge=0)) - mass.fast_mass('') + settings.getfloat('modifications', 'protein nterm cleavage') + settings.getfloat('modifications', 'protein cterm cleavage')
        for modif in self.modifications:
            self.pmass += modif['mass']
            if modif['position'] not in [0, len(self.sequence) + 1]:
                aminoacid = self.sequence[modif['position'] - 1]
                self.pmass -= mass.std_aa_mass[aminoacid]
            else:
                if modif['position'] == 0:
                    self.pmass -= settings.getfloat('modifications', 'protein nterm cleavage')
                else:
                    self.pmass -= settings.getfloat('modifications', 'protein cterm cleavage')

        self.mz = (mass_exp + pcharge * 1.007276) / pcharge
        self.mass_exp = mass_exp
        self.modified_peptide()
        self.RT_exp = RT_exp
        self.RT_predicted = False
        self.evalue = float(evalue)
        self.parentproteins = []
        self.massdiff = float(mass_exp) - float(self.pmass)
        self.num_missed_cleavages = dict()
        self.mc = mc
        self.note = note
        self.note2 = ''
        self.note3 = ''
        self.possible_mass = []
        self.protscore2 = 1
        self.peptscore = 1
        self.peptscore2 = 1
        self.spectrum = spectrum
        self.spectrum_mz = None
        self.fragment_mt = None
        self.sumI = sumI

    def theor_spectrum(self, types=('b', 'y'), maxcharge=None, **kwargs):
        peaks = {}
        maxcharge = max(self.pcharge, 1)
        for ion_type in types:
            ms = []
            for i in range(1, len(self.modified_sequence)):
                if self.modified_sequence[i - 1] in parser.std_amino_acids and self.modified_sequence[i] != '-':
                    for charge in range(1, maxcharge + 1):
                        if ion_type[0] in 'abc':
                            ms.append(mass.fast_mass2(
                                str(self.modified_sequence)[:i], ion_type=ion_type, charge=charge,
                                **kwargs))
                        else:
                            ms.append(mass.fast_mass2(
                                str(self.modified_sequence)[i:], ion_type=ion_type, charge=charge,
                                **kwargs))
            marr = np.array(ms)
            marr.sort()
            peaks[ion_type] = marr
        return peaks

    def get_missed_cleavages(self, protease='trypsin'):
        if protease not in self.num_missed_cleavages:
            self.num_missed_cleavages[protease] = sum(1 for x in re.finditer(protease, self.sequence))
        return self.num_missed_cleavages[protease]

    def count_modifications(self, label):
        if ']' in label:
            naa = 1
            nmods = self.modified_sequence.count('-' + label[:-1])
        elif '[' in label:
            naa = 1
            nmods = self.modified_sequence.count(label[:-1] + '-')
        else:
            naa = self.modified_sequence.count(label[-1])
            nmods = self.modified_sequence.count(label)
        if naa:
            return float(nmods) / naa
        else:
            return -1.0

    def get_median_fragment_mt(self, settings=None):
        if not self.fragment_mt and len(self.spectrum_mz) > 1:
            temp = settings.get('fragment mass', 'ion types')
            ion_types = (x.strip() for x in temp.split(','))
            acc = settings.getfloat('fragment mass', 'mass accuracy')
            spectrum_mz = copy(self.spectrum_mz)
            int_array = copy(self.spectrum_i)
            int_array = int_array / int_array.max() * 100
            i = int_array > int_array.max() / 100
            spectrum_mz = spectrum_mz[i]
            theor = self.theor_spectrum(types=ion_types, aa_mass=self.aa_mass)
            spectrum_KDTree = cKDTree(spectrum_mz.reshape((spectrum_mz.size, 1)))
            dist_total = np.array([])
            for fragments in theor.values():
                n = fragments.size
                dist, ind = spectrum_KDTree.query(fragments.reshape((n, 1)),
                    distance_upper_bound=acc)
                dist_total = np.append(dist_total, dist[dist != np.inf])
            if dist_total.size:
                self.fragment_mt = np.median(dist_total)
            else:
                self.fragment_mt = acc
        return self.fragment_mt

    def mass_diff(self):
        """Calculates a difference between theoretical and experimental masses. Takes into account an isotope mass difference error"""
        return (self.massdiff - round(self.massdiff, 0) * 1.0033548378) / (self.pmass - round(self.massdiff, 0) * 1.0033548378) * 1e6

    def modified_peptide(self):
        def add_modification(arg, term=None):
            i = ''
            done_flag = 0
            while 1:
                for x in lowercase:
                    if i + x not in self.modification_list.values():
                        self.modification_list[arg] = i + x
                        if term and term == 'c':
                             self.modification_list[arg] = '-' + self.modification_list[arg]
                        elif term and term == 'n':
                             self.modification_list[arg] += '-'
                        else:
                            print 'label for %s modification is missing in parameters, using %s label' % (arg, self.modification_list[arg])
                        done_flag = 1
                        break
                if not done_flag:
                    if i and i[-1] != lowercase[-1]:
                        i = i[:-1] + lowercase.index(i[-1] + 1)
                    else:
                        i += lowercase[0]
                else:
                    break

        def get_modification(arg):
            if arg.isdigit():
                if arg not in self.modification_list:
                    add_modification(arg)
                for modif in self.modifications:
                    if int(modif['mass']) == int(arg):
                        self.aa_mass[self.modification_list[arg] + self.sequence[modif['position'] - 1]] = float(modif['mass'])
                return self.modification_list[arg]
            else:
                return arg

        self.modified_sequence = ''
        stack = []
        tcode = re.split('\[|\]', self.modified_code)
        for elem in ''.join(map(get_modification, tcode))[::-1]:
            if elem.islower():
                stack.append(elem)
            else:
                self.modified_sequence += elem
                while stack:
                    self.modified_sequence += stack.pop(0)
        while stack:
            self.modified_sequence += stack.pop(0)
        self.modified_sequence = self.modified_sequence[::-1]
        for idx, elem in enumerate(tcode):
            flag_v = 0
            while idx - flag_v - 1 >= 0:
                if not tcode[idx - flag_v - 1]:
                    flag_v += 2
                else:
                    break
            if elem.isdigit() and not tcode[idx + 1] and idx != len(tcode) - 2 and tcode[idx + 2].isdigit():
                self.aa_mass[self.modification_list[elem] + self.modification_list[tcode[idx + 2]]] = self.aa_mass[self.modification_list[elem] + tcode[idx - 1 - flag_v][-1]] + self.aa_mass[self.modification_list[tcode[idx + 2]] + tcode[idx - 1 - flag_v][-1]] - mass.std_aa_mass[tcode[idx - 1 - flag_v][-1]]
        for modif in self.modifications:
            if modif['position'] == 0:
                try:
                    self.modified_sequence = self.modification_list[str(int(modif['mass']))] + self.modified_sequence
                except:
                    add_modification(str(int(modif['mass'])), term='n')
                    print 'label for %s nterm modification is missing in parameters, using %s label' % (str(int(modif['mass'])), self.modification_list[str(int(modif['mass']))])
                    self.aa_mass[self.modification_list[str(int(modif['mass']))]] = float(modif['mass'])
                    self.modified_sequence = self.modification_list[str(int(modif['mass']))] + self.modified_sequence
            elif modif['position'] == len(self.sequence) + 1:
                try:
                    self.modified_sequence = self.modified_sequence + self.modification_list[str(int(modif['mass']))]
                except:
                    add_modification(str(int(modif['mass'])), term='c')
                    print 'label for %s cterm modification is missing in parameters, using label %s' % (str(int(modif['mass'])), self.modification_list[str(int(modif['mass']))])
                    self.aa_mass[self.modification_list[str(int(modif['mass']))]] = float(modif['mass'])
                    self.modified_sequence = self.modified_sequence + self.modification_list[str(int(modif['mass']))]