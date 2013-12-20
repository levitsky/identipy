import sys
import re
import string
from string import punctuation
from pyteomics import parser
from pyteomics.parser import cleave, expasy_rules
from pyteomics import pepxml, electrochem, mass
from pyteomics.auxiliary import linear_regression
from pyteomics import achrom
import numpy as np
from scipy.stats import scoreatpercentile
from copy import copy
from scipy.spatial import cKDTree
try:
    from configparser import RawConfigParser
except ImportError:
    from ConfigParser import RawConfigParser

def get_aa_mass(settings):
    aa_mass = mass.std_aa_mass.copy()
#    print aa_mass
    fmods = settings.get('modifications', 'fixed')
    if fmods:
        for mod in re.split(r'[,;]\s*', fmods):
            m, aa = parser._split_label(mod)
            aa_mass[aa] += settings.getfloat('modifications', m)
    vmods = settings.get('modifications', 'variable')
    if vmods:
        mods = [parser._split_label(l) for l in re.split(r',\s*', vmods)]
        for (mod, aa), char in zip(mods, punctuation):
            aa_mass[char] = aa_mass[aa] + settings.getfloat('modifications', mod)
    aa_mass['|'] = 42.0106
    return aa_mass


def get_settings(fname=None, default_name='default.cfg'):
#os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'default.cfg')
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


def get_descriptors(dictonary):
    descriptors = []
    for descr in dictonary.values():
        descriptors.append(Descriptor(name=descr['name'], ))


class Descriptor():
    def __init__(self, name, formula, group='A', binsize='auto'):
        self.name = name
        self.formula = formula
        self.group = group
        self.binsize = binsize
        self.normK = None

    def get_array(self, peptides):
        code = """[%s for peptide in peptides.peptideslist]""" % (self.formula)
        return eval(code)

    def get_binsize(self, peptides=None):
        if self.binsize == 'auto':
            if self.group == 'B':
                return 0.25
            else:
                return FDbinSize(self.get_array(peptides))
        else:
            return self.binsize


class PeptideList:
    def __init__(self, settings=None):
        self.peptideslist = []
        self.calibrate_coeff = None
        self.pepxml_type = ''
        self.RC = False
        self.settings = settings

    def get_from_pepxmlfile(self, pepxmlfile, min_charge=1, max_charge=0, max_rank=1):
        for line in open(pepxmlfile, 'r').readlines():
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

        for record in pepxml.read(pepxmlfile):
            if 'search_hit' in record.keys():
                if int(min_charge) <= int(record['assumed_charge']) and (int(record['assumed_charge']) <= int(max_charge) or not max_charge):
                    for k in range(min(len(record['search_hit']), max_rank)):
                        sequence = record['search_hit'][k]['peptide']
                        if all(aminoacid not in ['B', 'X', 'J', 'Z', 'U', 'O', '.'] for aminoacid in sequence):
                            start_scan = record['start_scan']
                            modified_code = record['search_hit'][k]['modified_peptide']
                            modifications = record['search_hit'][k]['modifications']
                            prev_aa = record['search_hit'][k]['proteins'][0]['peptide_prev_aa']
                            next_aa = record['search_hit'][k]['proteins'][0]['peptide_next_aa']
                            try:
                                evalue = record['search_hit'][k]['search_score']['expect']
                            except:
                                try:
                                    evalue = 1.0 / float(record['search_hit'][k]['search_score']['ionscore'])
                                except IOError:
                                    'Cannot read e-value!'
                            try:
                                hyperscore = record['search_hit'][k]['search_score']['hyperscore']
                                nextscore = record['search_hit'][k]['search_score']['nextscore']
                            except:
                                hyperscore = 1
                                nextscore = 1
                            massdiff = record['search_hit'][k]['massdiff']
                            spectrum = record['spectrum']
                            rank = k + 1
                            pcharge = record['assumed_charge']
                            mass_exp = record['precursor_neutral_mass']

                            pept = Peptide(sequence=sequence, settings=self.settings, modified_code=modified_code, evalue=evalue, massdiff=massdiff, spectrum=spectrum, rank=rank, pcharge=pcharge, mass_exp=mass_exp, hyperscore=hyperscore, nextscore=nextscore, prev_aa=prev_aa, next_aa=next_aa, start_scan=start_scan, modifications=modifications)
                            try:
                                pept.RT_exp = float(record['retention_time_sec']) / 60
                            except:
                                try:
                                    pept.RT_exp = float(record['spectrum'].split(',')[2].split()[0])
                                except:
                                    pept.RT_exp = 0
                                    print 'RT_exp read error'

                            decoy_tags = [':reversed', 'DECOY_', 'rev_', 'Random sequence.']
                            if any([all([all([key in protein.keys() and isinstance(protein[key], str) and not protein[key].startswith(tag) and not protein[key].endswith(tag) for key in ['protein', 'protein_descr']]) for tag in decoy_tags]) for protein in record['search_hit'][k]['proteins']]):
                                pept.note = 'target'
                            else:
                                pept.note = 'decoy'

                            def get_dbname(prot, pepxml_type='tandem'):
                                if pepxml_type != 'omssa':
                                    try:
                                        if not any(prot['protein'].startswith(tag) for tag in ['sp', 'tr', 'DECOY_sp', 'DECOY_tr']):
                                            if any(prot['protein_descr'].startswith(tag) for tag in ['SWISS-PROT:', 'TREMBL:']):
                                                return prot['protein']
                                            if '|' not in prot['protein']:
                                                return prot['protein']+' '+prot['protein_descr']
                                            else:
                                                return prot['protein']+'>'+prot['protein_descr']
                                        else:
                                            return prot['protein'].split('|')[1]
                                    except:
                                        return prot['protein'].split('_')[0]
                                else:
                                    return prot['protein_descr'].split('|')[1]

                            for prot in record['search_hit'][k]['proteins']:
                                if get_dbname(prot, self.pepxml_type) not in [protein.dbname for protein in pept.parentproteins]:
                                    pept.parentproteins.append(Protein(dbname=get_dbname(prot, self.pepxml_type)))

                            if len(pept.parentproteins):
                                self.peptideslist.append(pept)

    def modified_peptides(self):
        for peptide in self.peptideslist:
            peptide.modified_peptide()

    def get_RC(self):
        seqs = [pept.modified_sequence for pept in self.peptideslist]
        RTexp = [pept.RT_exp for pept in self.peptideslist]
        RC_def = achrom.RCs_gilar_rp
        xdict = {}
        for key, val in RC_def['aa'].items():
            xdict[key] = [val, None]
        RC_dict = achrom.get_RCs_vary_lcp(seqs, RTexp)
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
        if 'C' not in RC_dict['aa'].keys():
            RC_dict['aa']['C'] = RC_dict['aa']['C*']
        self.RC = RC_dict

    def calc_RT(self, calibrate_coeff=(1, 0, 0, 0), RTtype='achrom'):
        if RTtype == 'ssrcalc':
            ps = list(set([peptide.sequence for peptide in self.peptideslist]))
            SSRCalc_RTs = SSRCalc.calculate_RH(ps[:], pore_size=100, ion_pairing_agent='FA')

        for peptide in self.peptideslist:
            if RTtype == 'achrom':
                peptide.RT_predicted = achrom.calculate_RT(peptide.modified_sequence, self.RC, raise_no_mod=False)
            elif RTtype == 'ssrcalc':
                SSRCalc_RT = SSRCalc_RTs[peptide.sequence]
                if SSRCalc_RT is not None:
                    peptide.RT_predicted = float(SSRCalc_RT) * calibrate_coeff[0] + calibrate_coeff[1]
                else:
                    peptide.RT_predicted = 0
                    print 'SSRCalc error'
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
            if peptide.sequence not in peptides_added.keys():
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

    def filter_unknown_aminoacids(self):
        j = len(self.peptideslist) - 1
        while j >= 0:
            if any(aminoacid in ['B', 'X', 'J', 'Z', 'U', 'O', '.']
                for aminoacid in self.peptideslist[j].sequence):
                    self.peptideslist.pop(j)
            j -= 1

    def filter_modifications(self):
        j = len(self.peptideslist) - 1
        while j >= 0:
            if (self.peptideslist[j].modified_code.count('[') - self.peptideslist[j].modified_code.count('[160]')) != 0:
                self.peptideslist.pop(j)
            j -= 1

    def filter_modifications_test(self):
        j = len(self.peptideslist) - 1
        while j >= 0:
            if self.peptideslist[j].modified_code.count('[') != 0:
                self.peptideslist.pop(j)
            j -= 1

    def filter_decoy(self):
        j = len(self.peptideslist) - 1
        while j >= 0:
            if self.peptideslist[j].note == 'decoy':
                self.peptideslist.pop(j)
            j -= 1

    def filter_titin(self):
        j = len(self.peptideslist) - 1
        while j >= 0:
            if any([prot.dbname == 'Q8WZ42' for prot in self.peptideslist[j].parentproteins]):
                self.peptideslist.pop(j)
            j -= 1

#    def calculate_proteins_score(self):
#        prots = {}
#        for peptide in self.peptideslist:
#            for prot in peptide.parentproteins:
#                if prot.dbname in prots.keys():
#                    prots[prot.dbname] += float(len(peptide.sequence))
#                else:
#                    prots[prot.dbname] = float(len(peptide.sequence))
#        for peptide in self.peptideslist:
#            for prot in peptide.parentproteins:
#                prot.score = float(prots[prot.dbname]) / len(prot.sequence)

    def filter_evalue_new(self, FDR=1, useMP=True, k=0):
        "A function for filtering PSMs by e-value and MP-score with some FDR"
        target_evalues, decoy_evalues = np.array([]), np.array([])
        for peptide in self.peptideslist:
            if peptide.note == 'target':
                target_evalues = np.append(target_evalues, float(peptide.evalue))
            elif peptide.note == 'decoy':
                decoy_evalues = np.append(decoy_evalues, float(peptide.evalue))
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

        best_cut_peptscore = 1
        if useMP:
            target_peptscores, decoy_peptscores = np.array([]), np.array([])
            for peptide in self.peptideslist:
                if peptide.evalue >= best_cut_evalue:
                    if peptide.note == 'target':
                        target_peptscores = np.append(target_peptscores, float(peptide.peptscore))
                    elif peptide.note == 'decoy':
                        decoy_peptscores = np.append(decoy_peptscores, float(peptide.peptscore))
            target_peptscores = np.sort(target_peptscores)[::-1]
            real_FDR = 0
            for cut_peptscore in target_peptscores:
                counter_target = target_peptscores[target_peptscores >= cut_peptscore].size
                counter_decoy = decoy_peptscores[decoy_peptscores >= cut_peptscore].size
                if counter_target and (float(counter_decoy) / float(counter_target)) * 100 <= float(FDR) / (1 + k):
                    best_cut_peptscore = cut_peptscore
                    real_FDR = round(float(counter_decoy) / float(counter_target) * 100, 1)
            print real_FDR, best_cut_peptscore, 'MP score'
        j = len(self.peptideslist) - 1
        while j >= 0:
            if not useMP or self.peptideslist[j].peptscore < best_cut_peptscore:
                if self.peptideslist[j].evalue > best_cut_evalue:
                    self.peptideslist.pop(j)
            j -= 1
        self.filter_decoy()
        return (best_cut_evalue, best_cut_peptscore)


class Protein:
    def __init__(self, dbname, pcharge=0, description='Unknown', sequence='Unknown', note=''):
        self.dbname = dbname
        self.pcharge = pcharge
        self.description = description
        self.PE = 0
        self.sequence = sequence
        self.peptides_exp = []
        self.peptides_theor = []
        self.pmass = 0
        self.note = note
        self.dbname2 = 'unknown'
        self.score = 0

    def get_mass(self):
        if self.sequence != 'Unknown':
            if any(aminoacid in ['B', 'X', 'J', 'Z', 'U', 'O', '.']
                for aminoacid in self.sequence):
                    self.pmass = 0
            else:
                self.pmass = float(mass.calculate_mass(sequence=self.sequence, charge=self.pcharge))

        else:
            print 'Unknown sequence'


class Peptide:
    def __init__(self, sequence, settings, modified_code='', pcharge=0, RT_exp=False, evalue=0, protein='Unkonwn', massdiff=0, note='unknown', spectrum='', rank=1, mass_exp=0, hyperscore=0, nextscore=0, prev_aa='X', next_aa='X', start_scan=0, modifications=[]):
        self.sequence = sequence
        self.modified_code = modified_code
        self.modified_sequence = sequence
        self.modifications = modifications
        self.pcharge = int(pcharge)
        self.nomodifications = 0
        self.unknown_mods_counter = 0 
        self.missed_modifications = []

        self.pmass = float(mass.calculate_mass(sequence=self.sequence, charge=0))
        for modif in self.modifications:
            self.pmass += modif['mass']
            if modif['position'] not in [0, len(self.sequence) + 1]:
                aminoacid = self.sequence[modif['position'] - 1]
                self.pmass -= mass.std_aa_mass[aminoacid]
            else:
                if modif['position'] == 0:
                    self.pmass -= 1.00782503207
                else:
                    self.pmass -= 17.002739651629998
        self.pI = electrochem.pI(self.sequence)

        for mod in self.modifications:
            if mod['position'] == 0:
#                self.modified_code = '[' + str(int(mod['mass'])) + ']' + self.modified_code
                self.sequence = '|' + self.sequence
                break

        self.mz = (mass_exp + pcharge * 1.007276) / pcharge
        self.mass_exp = mass_exp
        self.modified_peptide(settings)
        self.RT_exp = RT_exp
        self.RT_predicted = False
        self.evalue = float(evalue)
        self.start_scan = int(start_scan)
        self.parentprotein = protein
        self.parentproteins = []
        self.massdiff = float(mass_exp) - float(self.pmass)
        """
        if float(massdiff) != 0:
            self.massdiff = float(massdiff)
        else:
            if float(mass_exp) - float(self.pmass) < 0.01:
                self.massdiff = float(mass_exp) - float(self.pmass)
            else:
                self.massdiff = float(mass_exp) - float(self.pmass)

#                self.massdiff = 0
        """
        self.num_missed_cleavages = dict()
        self.note = note
        self.note2 = ''
        self.note3 = ''
        self.possible_mass = []
        self.protscore = 1
        self.protscore2 = 1
        self.peptscore = 1
        self.peptscore2 = 1
        self.spectrum = spectrum
        self.spectrum_mz = None
        self.fragment_mt = None
        self.rank = rank
        self.concentration = 1
        self.solubility = 0
        self.hyperscore = float(hyperscore)
        self.nextscore = float(nextscore)
        self.prev_aa = prev_aa
        self.next_aa = next_aa

    def theor_spectrum(self, types=('b', 'y'), maxcharge=None, **kwargs):
        peaks = {}
    #    maxcharge = 3
        maxcharge = self.pcharge
        for ion_type in types:
            ms = []
            for i in range(1, len(self.sequence) - 1):
                for charge in range(1, maxcharge + 1):
                    if ion_type[0] in 'abc':
                        ms.append(mass.fast_mass(
                            str(self.sequence)[:i], ion_type=ion_type, charge=charge,
                            **kwargs))
                    else:
                        ms.append(mass.fast_mass(
                            str(self.sequence)[i:], ion_type=ion_type, charge=charge,
                            **kwargs))
            marr = np.array(ms)
            marr.sort()
            peaks[ion_type] = marr
        return peaks

    def get_missed_cleavages(self, protease='trypsin'):
        if protease not in self.num_missed_cleavages.keys():
#            self.num_missed_cleavages = len(cleave(self.sequence, expasy_rules[protease], 0)) - 1
            self.num_missed_cleavages[protease] = len(cleave(self.sequence, protease, 0)) - 1
        return self.num_missed_cleavages[protease]

    def count_modifications(self, labels):
        nmods = 0
        for mod in labels:
#            print mod
            nmods += self.modified_sequence.count(mod)
        return nmods

    def count_modifications_new(self, label):
#        print label
        nmods = self.modified_sequence.count(label)
        naa = self.modified_sequence.count(parser._split_label(label)[1])
        if naa:
            return float(nmods) / naa
        else:
            return -1.0

    def get_fragment_mt(self, settings=None):
        if not self.fragment_mt:
            temp = settings.get('fragment mass', 'ion types')
            ion_types = (x.strip() for x in temp.split(','))
            acc = settings.getfloat('fragment mass', 'mass accuracy')
#            acc = fragment_mass_acc
#            acc = 0.5
            spectrum_mz = copy(self.spectrum_mz)
#            idx = np.nonzero(spectrum_mz >= 150)
#            spectrum_mz = spectrum_mz[idx]
            theor = self.theor_spectrum(types=ion_types, aa_mass=get_aa_mass(settings))
            spectrum_KDTree = cKDTree(spectrum_mz.reshape((spectrum_mz.size, 1)))

            dist_total = np.array([])
            for fragments in theor.values():
                n = fragments.size
                dist, ind = spectrum_KDTree.query(fragments.reshape((n, 1)),
                    distance_upper_bound=acc)
                mask = (dist != np.inf)
                dist_total = np.append(dist_total, dist[dist != np.inf])
            if dist_total.size:
                self.fragment_mt = np.median(dist_total)
#                self.fragment_mt = float(scoreatpercentile(dist_total, 75) - scoreatpercentile(dist_total, 25))
#                 self.fragment_mt = scoreatpercentile(dist_total, 80)
            else:
                self.fragment_mt = acc
        return self.fragment_mt

    def get_fragment_mt_ppm(self, settings=None):
        if not self.fragment_mt:
            temp = settings.get('fragment mass', 'ion types')
            ion_types = (x.strip() for x in temp.split(','))
            acc = settings.getfloat('fragment mass', 'mass accuracy')
#            acc = fragment_mass_acc
#            acc = 0.5
            spectrum_mz = copy(self.spectrum_mz)
#            idx = np.nonzero(spectrum_mz >= 150)
#            spectrum_mz = spectrum_mz[idx]
            theor = self.theor_spectrum(types=ion_types, aa_mass=get_aa_mass(settings))
            spectrum_KDTree = cKDTree(spectrum_mz.reshape((spectrum_mz.size, 1)))

            dist_total = np.array([])
            for fragments in theor.values():
                n = fragments.size
                dist, ind = spectrum_KDTree.query(fragments.reshape((n, 1)),
                    distance_upper_bound=acc)
                mask = (dist != np.inf)
#                temp = copy(dist[dist != np.inf])
                for j, v in enumerate(dist):
                    if v != np.inf:
#                        print j
    #                    print dist[j], '!'
                        dist[j] = float(dist[j]) * 1e6 / fragments[j]
    #                    print dist[j], '!!'
                dist_total = np.append(dist_total, dist[dist != np.inf])
            if dist_total.size:
                self.fragment_mt = np.median(dist_total)
            else:
                self.fragment_mt = acc
        return self.fragment_mt

    def mass_diff(self):
        """Calculates a difference between theoretical and experimental masses. Takes into account an isotope mass difference error"""
        return (self.massdiff - round(self.massdiff, 0) * 1.0033548378) / (self.pmass - round(self.massdiff, 0) * 1.0033548378) * 1e6

    def modified_peptide(self, settings):
        modifications = {}
        fmods = settings.get('modifications', 'fixed')
        if fmods:
            for mod in re.split(r'[,;]\s*', fmods):
                m, aa = parser._split_label(mod)
                modifications[str(int(mass.std_aa_mass[aa] + settings.getfloat('modifications', m)))] = m
        vmods = settings.get('modifications', 'variable')
        if vmods:
            mods = [parser._split_label(l) for l in re.split(r',\s*', vmods)]
            for (mod, aa), char in zip(mods, punctuation):
                modifications[str(int(mass.std_aa_mass[aa] + settings.getfloat('modifications', mod)))] = mod
        i = 0
        start = 0
        end = 0
        self.modified_sequence = str(self.modified_code)
        while i <= len(self.modified_sequence) - 1:
            if self.modified_sequence[i] == '[':
                start = i
            elif self.modified_sequence[i] == ']':
                end = i
                if self.modified_sequence[start + 1:end] not in modifications.keys():
                    self.missed_modifications.append(self.modified_sequence[start + 1:end])
#                    print 'missed modification with mass = %s' % (float(self.modified_sequence[start + 1:end]), )
                    mod_label = 'u'
#                    print modifications.values()
                    while mod_label in modifications.values():
                        mod_label += 'u'
                    mod_label = ''
                    temp_val = int(self.modified_sequence[start + 1:end])
                    maxV = len(string.lowercase)
                    while temp_val > maxV:
                        mod_label += string.lowercase[-1]
                        temp_val -= maxV
                    mod_label += string.lowercase[temp_val - 1]
                        
                    modifications[self.modified_sequence[start + 1:end]] = mod_label
                    for md_mass in self.modifications:
                        if int(md_mass['mass']) == int(self.modified_sequence[start + 1:end]):
                            temp_mass = md_mass['mass']
                            break
                    char = mod_label + self.modified_sequence[start-1]
#                    print mass.std_aa_mass
                    mass.std_aa_mass[char] = temp_mass
#                    print mass.std_aa_mass
                    modifications[int(temp_mass)] = mod_label
#                    print char, temp_mass, mod_label, '!@#!#@!$!@$!@$!@$!$@!$@!$@!@'
                    label = mod_label
                    #label = ''
                    self.unknown_mods_counter += 1
                else:
                    label = modifications[self.modified_sequence[start + 1:end]]
#                    print label, 'ok'
                self.modified_sequence = self.modified_sequence[:start - 1] + label + self.modified_sequence[start - 1] + self.modified_sequence[end + 1:]
                i = 0
            i += 1
