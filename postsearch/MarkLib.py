from pyteomics.parser import cleave, expasy_rules
from pyteomics import pepxml, electrochem, mass
from pyteomics.auxiliary import linear_regression
from pyteomics import achrom
import SSRCalc
import numpy as np
from freedman_diaconis import FDbinSize
from scipy.stats import scoreatpercentile

modifications = {160: 'cam', 147: 'ox', 181: 'p', 167: 'p', 243: 'p'}
constant_modifications = {'cam': ['C']}
potential_modifications = {'ox': ['M'], 'p': ['S', 'T', 'Y']}

mass.std_aa_comp['p'] = mass.Composition('HPO3')
mass.std_aa_comp['ox'] = mass.Composition('O')
mass.std_aa_comp['cam'] = mass.Composition('H3C2NO')


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

    def get_array(self, peptides):
        code = """[%s for peptide in peptides.peptideslist]""" % (self.formula)
        return eval(code)

    def get_binsize(self, peptides=None):
        if self.binsize == 'auto':
            if self.group == 'A':
                return FDbinSize(self.get_array(peptides))
            else:
                return 0.25
        else:
            return self.binsize


class PeptideList:
    def __init__(self):
        self.peptideslist = []
        self.calibrate_coeff = None
        self.pepxml_type = ''
        self.RC = False

    def get_from_txt(self, txtfile):
        for line in open(txtfile):
            sequence = line.split('\t')[0]
            modified_code = sequence
            evalue = float(line.split('\t')[1])
            massdiff = float(line.split('\t')[2])
            spectrum = 's'
            rank = 1
            pcharge = 100
            mass_exp = 10000000
            pept = Peptide(sequence=sequence, modified_code=modified_code, evalue=evalue, massdiff=massdiff, spectrum=spectrum, rank=rank, pcharge=pcharge, mass_exp=mass_exp)
            pept.RT_exp = float(line.split('\t')[3])
            pept.note = 'decoy'
            j = len(line.split('\t')) - 1
            while j > 3:
                p = line.split('\t')[j].strip()[1:-1]
                dbname = p.split(',')[0][1:-1]
                note = p.split(',')[1][2:-1]
                if note == 't':
                    pept.note = 'target'
                pept.parentproteins.append(Protein(dbname=dbname))
                j -= 1
            self.peptideslist.append(pept)

    def get_from_pepxmlfile(self, pepxmlfile, min_charge=1, max_charge=10, max_rank=1):
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
                if int(min_charge) <= int(record['assumed_charge']) <= int(max_charge):
                    for k in range(min(len(record['search_hit']), max_rank)):
                        sequence = record['search_hit'][k]['peptide']
                        if all(aminoacid not in ['B', 'X', 'J', 'Z', 'U', 'O'] for aminoacid in sequence):
                            modified_code = record['search_hit'][k]['modified_peptide']
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

                            pept = Peptide(sequence=sequence, modified_code=modified_code, evalue=evalue, massdiff=massdiff, spectrum=spectrum, rank=rank, pcharge=pcharge, mass_exp=mass_exp, hyperscore=hyperscore, nextscore=nextscore, prev_aa=prev_aa, next_aa=next_aa)
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
        seqs = [pept.sequence for pept in self.peptideslist]
        RTexp = [pept.RT_exp for pept in self.peptideslist]
        RC_def = achrom.RCs_yoshida
        xdict = {}
        for key, val in RC_def['aa'].items():
            xdict[key] = [val, None]
        RC_dict = achrom.get_RCs_vary_lcp(seqs, RTexp)
        for key, val in RC_dict['aa'].items():
            xdict[key][1] = val
        a, b, _, _ = linear_regression([x[0] for x in xdict.values() if x[1] != None], [x[1] for x in xdict.values() if x[1] != None])
        for key, x in xdict.items():
            if x[1] == None:
                x[1] = x[0] * a + b
            RC_dict['aa'][key] = x[1]
        self.RC = RC_dict

    def calc_RT(self, calibrate_coeff=(1, 0, 0, 0), RTtype='biolccc'):
        if RTtype == 'ssrcalc':
            ps = list(set([peptide.sequence for peptide in self.peptideslist]))
            SSRCalc_RTs = SSRCalc.calculate_RH(ps[:], pore_size=100, ion_pairing_agent='FA')

        for peptide in self.peptideslist:
            if RTtype == 'biolccc':
                peptide.RT_predicted = achrom.calculate_RT(peptide.sequence, self.RC)
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

    def get_calibrate_coeff(self, RTtype='biolccc'):
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
            if any(aminoacid in ['B', 'X', 'J', 'Z', 'U', 'O']
                for aminoacid in self.peptideslist[j].sequence):
                    self.peptideslist.pop(j)
            j -= 1

    def filter_modifications(self):
        j = len(self.peptideslist) - 1
        while j >= 0:
            if (self.peptideslist[j].modified_code.count('[') - self.peptideslist[j].modified_code.count('[160]')) != 0:
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

    def calculate_proteins_score(self):
        prots = {}
        for peptide in self.peptideslist:
            for prot in peptide.parentproteins:
                if prot.dbname in prots.keys():
                    prots[prot.dbname] += float(len(peptide.sequence))
                else:
                    prots[prot.dbname] = float(len(peptide.sequence))
        for peptide in self.peptideslist:
            for prot in peptide.parentproteins:
                prot.score = float(prots[prot.dbname]) / len(prot.sequence)

    def filter_evalue_new(self, FDR=1, useMP=True):
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
        print real_FDR, best_cut_evalue

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
                counter_target = target_peptscores[target_peptscores > cut_peptscore].size
                counter_decoy = decoy_peptscores[decoy_peptscores > cut_peptscore].size
                if counter_target and (float(counter_decoy) / float(counter_target)) * 100 <= float(FDR):
                    best_cut_peptscore = cut_peptscore
                    real_FDR = round(float(counter_decoy) / float(counter_target) * 100, 1)
            print real_FDR, best_cut_peptscore
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
            if any(aminoacid in ['B', 'X', 'J', 'Z', 'U']
                for aminoacid in self.sequence):
                    self.pmass = 0
            else:
                self.pmass = float(mass.calculate_mass(sequence=self.sequence, charge=self.pcharge))

        else:
            print 'Unknown sequence'


class Peptide:
    def __init__(self, sequence, modified_code='', pcharge=0, RT_exp=False, evalue=0, protein='Unkonwn', massdiff=0, note='unknown', spectrum='', rank=1, mass_exp=0, hyperscore=0, nextscore=0, prev_aa='X', next_aa='X'):
        self.sequence = sequence
        self.modified_code = modified_code
        self.modified_sequence = sequence
        self.pcharge = int(pcharge)
        self.nomodifications = 0
        self.pmass = float(mass.calculate_mass(sequence=self.sequence, charge=0))
        self.mz = (mass_exp + pcharge * 1.007276) / pcharge
#        self.modified_peptide()
        self.RT_exp = RT_exp
        self.RT_predicted = False
        self.evalue = float(evalue)
        self.parentprotein = protein
        self.parentproteins = []
        if float(massdiff) != 0:
            self.massdiff = float(massdiff)
        else:
            self.massdiff = float(mass_exp) - float(self.pmass)
        self.num_missed_cleavages = len(cleave(sequence, expasy_rules['trypsin'], 0)) - 1
        self.note = note
        self.note2 = ''
        self.possible_mass = []
        self.protscore = 1
        self.protscore2 = 1
        self.peptscore = 1
        self.peptscore2 = 1
        self.spectrum = spectrum
        self.rank = rank
        self.concentration = 1
        self.solubility = 0
        self.hyperscore = float(hyperscore)
        self.nextscore = float(nextscore)
        self.prev_aa = prev_aa
        self.next_aa = next_aa
        self.pI = electrochem.pI(self.sequence)

    def mass_diff(self):
        """Calculates a difference between theoretical and experimental masses. Takes into account an isotope mass difference error"""
        return (self.massdiff - round(self.massdiff, 0) * 1.0033548378) / (self.pmass - round(self.massdiff, 0) * 1.0033548378) * 1e6

    def modified_peptide(self):
        i = 0
        start = 0
        end = 0
        self.modified_sequence = str(self.modified_code)
        while i <= len(self.modified_sequence) - 1:
            if self.modified_sequence[i] == '[':
                start = i
            elif self.modified_sequence[i] == ']':
                end = i
                if float(self.modified_sequence[start + 1:end]) not in modifications.keys():
                    label = ''
                else:
                    label = modifications[float(self.modified_sequence[start + 1:end])]
                self.modified_sequence = self.modified_sequence[:start - 1] + label + self.modified_sequence[start - 1] + self.modified_sequence[end + 1:]
                i = 0
            i += 1
