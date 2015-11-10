import re
from pyteomics import mass, electrochem as ec, auxiliary as aux, parser, fasta
import sys
from itertools import combinations
import numpy as np
from multiprocessing import Queue, Process, cpu_count
from string import punctuation
from copy import copy
from ConfigParser import RawConfigParser

try:
    from pyteomics import cmass
except ImportError:
    cmass = mass

def normalize_mods(sequence, settings):
    leg = settings.get('misc', 'legend')
    if leg:
        for char in punctuation:
            if char in leg:
                sequence = sequence.replace(char, ''.join(leg[char]))
    return sequence

def custom_isoforms(peptide, variable_mods, maxmods=2):
    if not variable_mods:
        yield peptide
    else:
        to_char = variable_mods[-1][0]
        from_char = variable_mods[-1][1]
        sites = [s[0] for s in enumerate(peptide) if s[1] == from_char]
        for m in range(maxmods+1):
            for comb in combinations(sites, m):
                v = ''
                cc_prev = 0
                for cc in comb:
                    v = v + peptide[cc_prev:cc] + to_char
                    cc_prev = cc + 1
                v = v + peptide[cc_prev:]
                for z in custom_isoforms(v, variable_mods[:-1], maxmods=maxmods - m):
                    yield z

def preprocess_spectrum(spectrum, settings):
    spectrum = copy(spectrum)
    maxpeaks = settings.getint('scoring', 'maximum peaks')
    minpeaks = settings.getint('scoring', 'minimum peaks')
    dynrange = settings.getfloat('scoring', 'dynamic range')

    idx = np.nonzero(spectrum['m/z array'] >= settings.getfloat('search', 'product minimum m/z'))
    spectrum['intensity array'] = spectrum['intensity array'][idx]
    spectrum['m/z array'] = spectrum['m/z array'][idx]

    if minpeaks and spectrum['intensity array'].size < minpeaks:
        return None

    if dynrange:
        i = spectrum['intensity array'] > spectrum['intensity array'].max(
                ) / dynrange
        spectrum['intensity array'] = spectrum['intensity array'][i]
        spectrum['m/z array'] = spectrum['m/z array'][i]

    if maxpeaks and minpeaks > maxpeaks:
        raise ValueError('minpeaks > maxpeaks: {} and {}'.format(
            minpeaks, maxpeaks))
    if maxpeaks and spectrum['intensity array'].size > maxpeaks:
        i = np.argsort(spectrum['intensity array'])[-maxpeaks:]
        j = np.argsort(spectrum['m/z array'][i])
        spectrum['intensity array'] = spectrum['intensity array'][i][j]
        spectrum['m/z array'] = spectrum['m/z array'][i][j]
    
    if minpeaks and spectrum['intensity array'].size < minpeaks:
        return None
    return spectrum

def relative(unit):
    if unit == 'ppm':
        return True
    elif unit in {'Th', 'Da', 'amu'}:
        return False
    else:
        raise ValueError('Unrecognized precursor accuracy unit: ' + unit)

def set_mod_dict(settings):
    mods = settings.get('modifications', 'variable')
    if isinstance(mods, basestring):
        mods = mods.strip()
        mod_dict = {}
        legend = {}
        
        if mods:
            mods = [parser._split_label(l) for l in re.split(r',\s*', mods)]
            mods.sort(key=lambda x: len(x[0]), reverse=True)
            for mod, char in zip(mods, punctuation):
                legend[''.join(mod)] = char
                legend[char] = mod
            assert all(len(m) == 2 for m in mods), 'unmodified residue given'
            for mod, aa in mods:
                mod_dict.setdefault(mod, []).append(aa)
            settings.set('modifications', 'variable', mod_dict)
        settings.set('misc', 'legend', legend)

def get_enzyme(enzyme):
    if enzyme in parser.expasy_rules:
        return parser.expasy_rules.get(enzyme)
    else:
        try:
            enzyme = convert_tandem_cleave_rule_to_regexp(enzyme)
            return enzyme
        except:
            return enzyme

def convert_tandem_cleave_rule_to_regexp(cleavage_rule):

    def get_sense(c_term_rule, n_term_rule):
        if '{' in c_term_rule:
            return 'N'
        elif '{' in n_term_rule:
            return 'C'
        else:
            if len(c_term_rule) <= len(n_term_rule):
                return 'C'
            else:
                return 'N'

    def get_cut(cut, no_cut):
        aminoacids = set(parser.std_amino_acids)
        cut = ''.join(aminoacids & set(cut))
        if '{' in no_cut:
            no_cut = ''.join(aminoacids & set(no_cut))
            return cut, no_cut
        else:
            no_cut = ''.join(set(parser.std_amino_acids) - set(no_cut))
            return cut, no_cut

    out_rules = []
    for protease in cleavage_rule.split(','):
        protease = protease.replace('X', ''.join(parser.std_amino_acids))
        c_term_rule, n_term_rule = protease.split('|')
        sense = get_sense(c_term_rule, n_term_rule)
        if sense == 'C':
            cut, no_cut = get_cut(c_term_rule, n_term_rule)
        else:
            cut, no_cut = get_cut(n_term_rule, c_term_rule)

        if no_cut:
            if sense == 'C':
                out_rules.append('([%s](?=[^%s]))' % (cut, no_cut))
            else:
                out_rules.append('([^%s](?=[%s]))' % (no_cut, cut))
        else:
            if sense == 'C':
                out_rules.append('([%s])' % (cut, ))
            else:
                out_rules.append('(?=[%s])' % (cut, ))
    return '|'.join(out_rules)

class CustomRawConfigParser(RawConfigParser):
    def get(self, section, option):
        val = RawConfigParser.get(self, section, option)
        if isinstance(val, basestring):
            return val[::-1].split('|', 1)[-1][::-1]
        return val

    def get_choices(self, section, option):
        val = RawConfigParser.get(self, section, option)
        if isinstance(val, basestring) and len(val.split('|')) > 1:
            return val[::-1].split('|', 1)[0][::-1]
        else:
            return ''


def find_nearest(array, value):
    return (np.abs(np.array(array) - value)).argmin()


def get_info(spectrum, result, settings, aa_mass=None):
    'Returns neutral mass, charge state and retention time of the top candidate'
    if not aa_mass:
        aa_mass = get_aa_mass(settings)
    if 'params' in spectrum:
        RT = spectrum['params'].get('rtinseconds')
    else:
        RT = spectrum['scanList']['scan'][0]['scan start time']
    masses, states = zip(*neutral_masses(spectrum, settings))
    idx = find_nearest(masses, cmass.fast_mass(str(result['candidates'][0][1]), aa_mass=aa_mass))
    return (masses[idx], states[idx], RT)

def theor_spectrum(peptide, types=('b', 'y'), maxcharge=None, **kwargs):
    peaks = {}
    if not maxcharge:
        maxcharge = 1 + int(ec.charge(peptide, pH=2))
    for ion_type in types:
        ms = []
        for i in range(1, len(peptide) - 1):
            for charge in range(1, maxcharge + 1):
                if ion_type[0] in 'abc':
                    ms.append(cmass.fast_mass(
                        str(peptide)[:i], ion_type=ion_type, charge=charge,
                        **kwargs))
                else:
                    ms.append(cmass.fast_mass(
                        str(peptide)[i:], ion_type=ion_type, charge=charge,
                        **kwargs))
        marr = np.array(ms)
        marr.sort()
        peaks[ion_type] = marr
    return peaks


def neutral_masses(spectrum, settings):
    if 'params' in spectrum:
        exp_mass = spectrum['params']['pepmass'][0]
        charge = spectrum['params'].get('charge')
    else:
        ion = spectrum['precursorList']['precursor'][
                0]['selectedIonList']['selectedIon'][0]
        charge = ion.get('charge state')
        if charge is not None: charge = [int(charge)]
        exp_mass = ion['selected ion m/z']
    if isinstance(charge, str):
        states = [s for s in aux._parse_charge(charge, True)
                if settings.getint('search', 'minimum charge'
                    ) <= s <= settings.getint('search', 'maximum charge')]
    elif charge is None:
        states = range(settings.getint('search', 'minimum charge'),
            1 + settings.getint('search', 'maximum charge'))
    else:
        states = [c for c in charge if settings.getint('search', 'minimum charge'
                    ) <= c <= settings.getint('search', 'maximum charge')]

    states.sort()
    return zip((c * (exp_mass - mass.nist_mass['H+'][0][0])
            for c in states), states)


@aux.memoize(10)
def import_(name):
    """Import a function by name: module.function or
    module.submodule.function, etc. By default trying to find
    function name in identipy.scoring module.
    Return the function object."""

    print 'Trying to import', name
    try:
        mod, f = name.rsplit('.', 1)
        return getattr(__import__(mod, fromlist=[f]), f)
    except Exception as e:
        print(e)
        return getattr(__import__('identipy.scoring', fromlist=[name]), name)

def get_aa_mass(settings):
    if settings.has_option('misc', 'aa_mass'):
        return settings.get('misc', 'aa_mass')
    aa_mass = mass.std_aa_mass.copy()
    for k, v in settings.items('modifications'):
        if k not in {'fixed', 'variable'}:
            aa_mass[k] = float(v)
    fmods = settings.get('modifications', 'fixed')
    if fmods:
        for mod in re.split(r'[,;]\s*', fmods):
            m, aa = parser._split_label(mod)
            aa_mass[aa] += settings.getfloat('modifications', m)
    vmods = settings.get('modifications', 'variable')
    if vmods:
        leg = settings.get('misc', 'legend')
        for p in punctuation:
            if p in leg:
                mod, aa = leg[p]
                aa_mass[p] = aa_mass[mod] + aa_mass[aa]
                aa_mass[mod+aa] = aa_mass[mod] + aa_mass[aa]
    return aa_mass

def cand_charge(result):
    mz = result['spectrum']['params']['pepmass'][0]
    m = np.array(map(lambda x: cmass.fast_mass(str(x)), result['candidates']['seq']))
    return np.round(m/mz).astype(np.int8)

def multimap(n, func, it, **kw):
    if n == 0:
        try:
            n = cpu_count()
        except NotImplementedError:
            n = 1
    if n == 1:
        for s in it:
            yield func(s, **kw)
    else:
        def worker(qin, qout):
            for spectrum in iter(qin.get, None):
                result = func(spectrum, **kw)
                qout.put(result)
        qin = Queue()
        qout = Queue()
        count = 0
        while True:
            procs = []
            for _ in range(n):
                p = Process(target=worker, args=(qin, qout))
                p.start()
                procs.append(p)
            for s in it:
                qin.put(s)
                count += 1
                if count > 500000:
                    print 'Loaded 500000 items. Ending cycle.'
                    break
            for _ in range(n):
                qin.put(None)

            if not count:
                print 'No items left. Exiting.'
                break

            while count:
                yield qout.get()
                count -= 1

            for p in procs:
                p.join()

            print 'Cycle finished.'

def allow_all(*args):
    return True


def get_RT(spectrum):
    if 'params' in spectrum:
        try:
            return float(spectrum['params']['rtinseconds']) / 60
        except:
            try:
                return float(spectrum['params']['title'].split(',')[-1].strip().split()[0])
            except:
                return 0
    return spectrum['scanList']['scan'][0]['scan start time']

def get_title(spectrum):
    if 'params' in spectrum:
        return spectrum['params']['title']

def get_output(results, settings):
    show_empty = settings.getboolean('output', 'show empty')
    score_threshold = settings.getfloat('output', 'score threshold')
    min_matched = settings.getint('output', 'minimum matched')
    num_candidates = settings.getint('output', 'candidates') or None
    acc_l = settings.getfloat('output', 'precursor accuracy left')
    acc_r = settings.getfloat('output', 'precursor accuracy right')
    rel = settings.get('output', 'precursor accuracy unit') == 'ppm'
    key = ['Th', 'ppm'][rel]
    count = 0
    for result in results:
        count += 1
        c = result['candidates']
        c = c[c['score'] > score_threshold]
        if min_matched:
            mask = np.array([
                c_[4]['match'] is not None and
                sum(m.sum() for m in c_[4]['match'].values()) >= min_matched
                for c_ in c], dtype=bool)
            c = c[mask]
        mask = np.array([-acc_l < c_[4]['mzdiff'][key] < acc_r for c_ in c], dtype=bool)
        c = c[mask]
        if (not c.size) and not show_empty:
            continue
        result['candidates'] = c[:num_candidates]
        yield result
    print 'Unfiltered results:', count


def write_pepxml(inputfile, settings, results):
    from lxml import etree
    from time import strftime
    from os import path

    if settings.has_option('output', 'path'):
        outpath = settings.get('output', 'path')
    else:
        outpath = path.dirname(inputfile)

    filename = path.join(outpath, path.splitext(path.basename(inputfile))[0] + path.extsep + 'pep' + path.extsep + 'xml')
    enzyme = settings.get('search', 'enzyme')
    search_engine = 'IdentiPy'
    database = settings.get('input', 'database')
    missed_cleavages = settings.getint('search', 'number of missed cleavages')
    fmods = settings.get('modifications', 'fixed')

    output = open(filename, 'w')
    line1 = '<?xml version="1.0" encoding="UTF-8"?>\n\
    <?xml-stylesheet type="text/xsl" href="pepXML_std.xsl"?>\n'
    output.write(line1)

    root = etree.Element('msms_pipeline_analysis')
    root.set("date", strftime("%Y:%m:%d:%H:%M:%S"))
    root.set("summary_xml", '')
    root.set("xmlns", 'http://regis-web.systemsbiology.net/pepXML')
    # TODO
    #root.set("xmlns:xsi", 'http://www.w3.org/2001/XMLSchema-instance')
    #root.set("xsi:schemaLocation", 'http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v117.xsd')

    child1 = etree.Element('msms_run_summary')
    child1.set("base_name", filename)
    child1.set("search_engine", search_engine)
    child1.set("raw_data_type", "raw")  # ?
    child1.set("raw_data", ".?")  # ?
    root.append(child1)

    child2 = etree.Element('sample_enzyme')
    child2.set('name', enzyme)
    child1.append(child2)

    child3 = etree.Element('specificity')
    child3.set("cut", "KR")
    child3.set("no_cut", "P")
    child3.set("sence", "C")

    child2.append(child3)

    child4 = etree.Element('search_summary')
    child4.set('base_name', filename)
    child4.set('search_engine', search_engine)
    child4.set('precursor_mass_type', 'monoisotopic')
    child4.set('fragment_mass_type', 'monoisotopic')
    child4.set('search_id', '1')

    child1.append(child4)

    child5 = etree.Element('search_database')
    child5.set('local_path', database)
    child5.set('type', 'AA')

    child4.append(copy(child5))

    child5 = etree.Element('enzymatic_search_constraint')
    child5.set('enzyme', enzyme)
    child5.set('max_num_internal_cleavages', str(missed_cleavages))
    child5.set('min_number_termini', '2')

    child4.append(copy(child5))

    results = [x for x in results if x['candidates'].size]
    results = list(get_output(results, settings))
    print 'Accumulated results:', len(results)
    pept_prot = dict()
    prots = dict()
    peptides = set()
    for x in results:
        peptides.update(re.sub(r'[^A-Z]', '', normalize_mods(x['candidates'][i][1], settings)) for i in range(
                settings.getint('output', 'candidates') or len(x['candidates'])))
    if settings.getboolean('input', 'add decoy'):
        decoy_method = settings.get('input', 'decoy method')
        decoy_prefix = settings.get('input', 'decoy prefix')
        f = fasta.decoy_db(database, mode=decoy_method, prefix=decoy_prefix)
    else:
        f = fasta.read(database)
    for desc, prot in f:
        dbinfo = desc.split(' ')[0]
        prots[dbinfo] = desc
        for pep in parser.cleave(prot, get_enzyme(enzyme), missed_cleavages):
            if pep in peptides:
                pept_prot.setdefault(pep, []).append(dbinfo)
    f.close()

    if settings.has_option('misc', 'aa_mass'):
        aa_mass = settings.get('misc', 'aa_mass')
    else:
        aa_mass = get_aa_mass(settings)

    vmods = set()
    variablemods =  settings.get('modifications', 'variable')
    if variablemods:
        for k, v in variablemods.items():
            for aa in v:
                vmods.add(k + aa)

    leg = {}
    if settings.has_option('misc', 'legend'):
        leg = settings.get('misc', 'legend')

    for idx, result in enumerate(results):
        if result['candidates'].size:
            tmp = etree.Element('spectrum_query')
            spectrum = result['spectrum']
            try:
                tmp.set('spectrum', spectrum['params']['title'])
            except:
                tmp.set('spectrum', spectrum['spectrum title'])
            tmp.set('start_scan', str(idx))  # ???
            tmp.set('end_scan', str(idx))  # ???
            tmp.set('index', str(idx))  # ???

            neutral_mass, charge_state, RT = get_info(spectrum, result, settings, aa_mass)
            tmp.set('precursor_neutral_mass', str(neutral_mass))
            tmp.set('assumed_charge', str(int(charge_state)))
            if RT:
                tmp.set('retention_time_sec', str(RT))

            tmp2 = etree.Element('search_result')
            result['candidates'] = result['candidates'][:len(result['e-values'])]

            flag = 1
            for i, candidate in enumerate(result['candidates']):
                if candidate[4]['match'] is None: break
                tmp3 = etree.Element('search_hit')
                tmp3.set('hit_rank', str(i + 1))
                mod_sequence = str(candidate[1])
                mod_sequence = normalize_mods(mod_sequence, settings)

                sequence = re.sub(r'[^A-Z]', '', mod_sequence)
                if sequence not in pept_prot:
                    flag = 0
                    break
                else:
                    tmp3.set('peptide', sequence)
                    tmp3.set('peptide_prev_aa', 'K')  # ???
                    tmp3.set('peptide_next_aa', 'K')  # ???
                    proteins = pept_prot[re.sub(r'[^A-Z]', '', sequence)]

                    tmp3.set('protein', prots[proteins[0]].split(' ', 1)[0])
                    tmp3.set('protein_descr', prots[proteins[0]].split(' ', 1)[1])

                    num_tot_proteins = len(proteins)
                    tmp3.set('num_tot_proteins', str(num_tot_proteins))
                    tmp3.set('num_matched_ions', '7')  # ???
                    tmp3.set('tot_num_ions', '7')  # ???
                    neutral_mass_theor = cmass.fast_mass(sequence, aa_mass=aa_mass)
                    tmp3.set('calc_neutral_pep_mass', str(neutral_mass_theor))
                    tmp3.set('massdiff', str(candidate[4]['mzdiff']['Th']))
                    tmp3.set('num_tol_term', '2')  # ???
                    tmp3.set('num_missed_cleavages', str(len(parser.cleave(sequence, get_enzyme(enzyme), 0)) - 1))
                    tmp3.set('is_rejected', '0')  # ???

                    if num_tot_proteins > 1:
                        for idx in range(len(proteins)):
                            if idx != 0:
                                tmp4 = etree.Element('alternative_protein')
                                tmp4.set('protein', prots[proteins[idx]].split(' ', 1)[0])
                                tmp4.set('protein_descr', prots[proteins[idx]].split(' ', 1)[1])
                                tmp3.append(copy(tmp4))

                    tmp4 = etree.Element('modification_info')
                    for idx, aminoacid in enumerate(parser.parse(mod_sequence)):
                        if aminoacid in fmods or aminoacid in vmods:
                            tmp5 = etree.Element('mod_aminoacid_mass')
                            tmp5.set('position', str(idx + 1))
                            tmp5.set('mass', str(aa_mass.get(aminoacid)))
                            tmp4.append(copy(tmp5))
                    tmp3.append(copy(tmp4))

                    tmp4 = etree.Element('search_score')
                    tmp4.set('name', 'hyperscore')
                    tmp4.set('value', str(candidate[0]))
                    tmp3.append(copy(tmp4))

                    tmp4 = etree.Element('search_score')
                    tmp4.set('name', 'nextscore')
                    tmp4.set('value', str(candidate[0]))
                    tmp3.append(copy(tmp4))

                    tmp4 = etree.Element('search_score')
                    tmp4.set('name', 'expect')
                    tmp4.set('value', str(result['e-values'][i]))
                    tmp3.append(copy(tmp4))

                    tmp4 = etree.Element('search_score')
                    tmp4.set('name', 'sumI')
                    tmp4.set('value', str(candidate[5]))
                    tmp3.append(copy(tmp4))

                    tmp2.append(copy(tmp3))
            if flag:
                tmp.append(copy(tmp2))
                child1.append(copy(tmp))

    s = etree.tostring(root, pretty_print=True)
    output.write(s)

    output.close()
