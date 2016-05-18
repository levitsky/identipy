import re
from pyteomics import mass, electrochem as ec, auxiliary as aux, fasta, mgf, mzml, parser
import sys
from itertools import combinations
import numpy as np
from multiprocessing import Queue, Process, cpu_count
from string import punctuation
from copy import copy
from ConfigParser import RawConfigParser
import tempfile
try:
    from pyteomics import cmass
except ImportError:
    cmass = mass

def custom_split_label(mod):
    j = 0
    while mod[j].islower():
        j += 1
    if j == 0:
        return mod[1:], '-', ']'
    if len(mod[j:]) > 1 and '[' in mod:
        return mod[:j], mod[j:].replace('[', ''), '['
    elif len(mod[j:]) > 1 and ']' in mod:
        return mod[:j], mod[j:].replace(']', ''), ']'
    elif len(mod[j:]) == 1:
        if mod.startswith('-'):
            return mod[:j], '-', ']'
        elif mod.endswith('-'):
            return mod[:j], '-', '['
        else:
            return mod[:j], mod[j:], ''

def iterate_spectra(fname):
    ftype = fname.rsplit('.', 1)[-1].lower()
    if ftype == 'mgf':
        with mgf.read(fname, read_charges=False) as f:
            for x in f:
                yield x
    elif ftype == 'mzml':
        with mzml.read(fname) as f:
            for x in f:
                if x['ms level'] > 1:
                    yield x
    else:
        raise ValueError('Unrecognized file type: {}'.format(ftype))

def peptide_gen(settings):
    prefix = settings.get('input', 'decoy prefix')
    enzyme = get_enzyme(settings.get('search', 'enzyme'))
    mc = settings.getint('search', 'number of missed cleavages')
    minlen = settings.getint('search', 'peptide minimum length')
    maxlen = settings.getint('search', 'peptide maximum length')
    for prot in prot_gen(settings):
        for pep in prot_peptides(prot[1], enzyme, mc, minlen, maxlen, is_decoy=prot[0].startswith(prefix)):
            yield pep

def prot_gen(settings):
    db = settings.get('input', 'database')
    add_decoy = settings.getboolean('input', 'add decoy')
    prefix = settings.get('input', 'decoy prefix')
    mode = settings.get('input', 'decoy method')

    read = [fasta.read, lambda f: fasta.decoy_db(f, mode=mode, prefix=prefix)][add_decoy]
    with read(db) as f:
        for p in f:
            yield p

seen_target = set()
seen_decoy = set()
def prot_peptides(prot_seq, enzyme, mc, minlen, maxlen, is_decoy):

    dont_use_fast_valid = parser.fast_valid(prot_seq)
    peptides = parser.cleave(prot_seq, enzyme, mc)
    for pep in peptides:
        plen = len(pep)
        if minlen <= plen <= maxlen + 2:
            forms = []
            if dont_use_fast_valid or pep in seen_target or pep in seen_decoy or parser.fast_valid(pep):
                # plen = len(pep)
                # if minlen <= plen <= maxlen:
                if plen <= maxlen:
                    forms.append(pep)
                if prot_seq[0] == 'M' and prot_seq.startswith(pep):
                    if minlen <= plen - 1 <= maxlen:
                        forms.append(pep[1:])
                    if minlen <= plen - 2 <= maxlen:
                        forms.append(pep[2:])
            for f in forms:
                if f not in seen_target and f not in seen_decoy:
                    if is_decoy:
                        seen_decoy.add(f)
                    else:
                        seen_target.add(f)
                    yield f


def normalize_mods(sequence, settings):
    leg = settings.get('misc', 'legend')
    if leg:
        for char in punctuation:
            if char in leg:
                if leg[char][2] == ']' and leg[char][1] == '-':
                    sequence = sequence.replace(char, '-' + leg[char][0])
                else:
                    sequence = sequence.replace(char, ''.join(leg[char][:2]))
    return sequence

def custom_isoforms(peptide, variable_mods, maxmods=2, nterm=False, cterm=False):
    if not variable_mods:
        yield peptide
    else:
        to_char = variable_mods[-1][0]
        from_char = variable_mods[-1][1]
        term = variable_mods[-1][2]
        sites = [s[0] for s in enumerate(peptide) if (from_char == '-' or s[1] == from_char) and (not term or (term == '[' and s[0] == 0) or (term == ']' and s[0] == len(peptide)-1))]
        for m in range(maxmods+1):
            for comb in combinations(sites, m):
                flag = 0
                flag2 = 0
                tmpnterm = True if nterm else False
                tmpcterm = True if cterm else False
                v = ''
                cc_prev = 0
                for cc in comb:
                    if from_char == '-':
                        if term == '[' and not nterm:
                            flag2 = 1
                            v += to_char
                            tmpnterm = True
                        elif term == ']' and not cterm:
                            v = v + peptide[cc_prev:cc+1] + to_char
                            tmpcterm = True
                        else:
                            flag = 1
                    else:
                        v = v + peptide[cc_prev:cc] + to_char
                    if not flag2:
                        cc_prev = cc + 1
                if not flag:
                    v = v + peptide[cc_prev:]
                    for z in custom_isoforms(v, variable_mods[:-1], maxmods=maxmods - m, nterm=tmpnterm, cterm=tmpcterm):
                        yield z

def preprocess_spectrum(spectrum, settings):
    spectrum = copy(spectrum)
    maxpeaks = settings.getint('scoring', 'maximum peaks')
    minpeaks = settings.getint('scoring', 'minimum peaks')
    dynrange = settings.getfloat('scoring', 'dynamic range')

    _, states = get_expmass(spectrum, settings)
    if not states:
        return None

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
            mods = [custom_split_label(l) for l in re.split(r',\s*', mods)]
            mods.sort(key=lambda x: len(x[0]), reverse=True)
            for mod, char in zip(mods, punctuation):
                legend[''.join(mod)] = char
                legend[char] = mod
            assert all(len(m) == 3 for m in mods), 'unmodified residue given'
            for mod, aa, term in mods:
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
    RT = get_RT(spectrum)
    masses, states = zip(*neutral_masses(spectrum, settings))
    idx = find_nearest(masses, cmass.fast_mass(str(result['candidates'][0][1]), aa_mass=aa_mass))
    return (masses[idx], states[idx], RT)

def theor_spectrum(peptide, types=('b', 'y'), maxcharge=None, reshape=False, **kwargs):
    peaks = {}
    pl = len(peptide) - 1
    if not maxcharge:
        maxcharge = 1 + int(ec.charge(peptide, pH=2))
    for charge in range(1, maxcharge + 1):
        for ion_type in types:
            ion_type_check = ion_type[0] in 'abc'
            ms = []
            for i in range(1, pl):
                if ion_type_check:
                    ms.append(cmass.fast_mass(
                        str(peptide)[:i], ion_type=ion_type, charge=charge,
                        **kwargs))
                else:
                    ms.append(cmass.fast_mass(
                        str(peptide)[i:], ion_type=ion_type, charge=charge,
                        **kwargs))
            marr = np.array(ms)
            if not reshape:
                marr.sort()
            else:
                n = marr.size
                marr = marr.reshape((n, 1))
            peaks[ion_type, charge] = marr
    return peaks

def get_expmass(spectrum, settings):
    maxcharge = settings.getint('search', 'maximum charge') or None
    mincharge = settings.getint('search', 'minimum charge') or None

    minu = settings.getint('search', 'minimum unknown charge')
    min_ucharge = minu if minu else mincharge
    maxu = settings.getint('search', 'maximum unknown charge')
    max_ucharge = maxu if maxu else maxcharge


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
                if (mincharge is None or s >= mincharge) and (maxcharge is None or s <= maxcharge)]
    elif charge is None:
        states = range(min_ucharge, 1 + max_ucharge)
    else:
        states = [c for c in charge if
            (mincharge is None or c >= mincharge) and (maxcharge is None or c <= maxcharge)]
    states.sort()
    return exp_mass, states


def neutral_masses(spectrum, settings):
    exp_mass, states = get_expmass(spectrum, settings)
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
    aa_mass['-'] = 0.0
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
                mod, aa, term = leg[p]
                if term == ']' and aa == '-':
                    aa_mass[p] = aa_mass[mod] + aa_mass[aa]
                    aa_mass[aa+mod] = aa_mass[mod] + aa_mass[aa]
                else:
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
            for item in iter(qin.get, None):
                result = func(item, **kw)
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
                if count > 5000000:
                    print 'Loaded 5000000 items. Ending cycle.'
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
                try:
                    return 60 * np.average([float(x) for x in spectrum['params']['title'].split('lution from: ')[-1].split(' period:')[0].split(' to ')])
                except:
                    return 0
    return spectrum['scanList']['scan'][0]['scan start time']

def get_title(spectrum):
    if 'params' in spectrum:
        return spectrum['params']['title']

def get_precursor_mz(spectrum):
    return spectrum['params']['pepmass'][0]

def get_output(results, settings):
    show_empty = settings.getboolean('output', 'show empty')
    score_threshold = settings.getfloat('output', 'score threshold')
    min_matched = settings.getint('output', 'minimum matched')
    num_candidates = settings.getint('output', 'candidates') or None
    try:
        acc_l = settings.getfloat('output', 'precursor accuracy left')
        acc_r = settings.getfloat('output', 'precursor accuracy right')
        rel = settings.get('output', 'precursor accuracy unit') == 'ppm'
    except ValueError:
        print 'Using [search] parameters for [output]'
        acc_l = settings.getfloat('search', 'precursor accuracy left')
        acc_r = settings.getfloat('search', 'precursor accuracy right')
        rel = settings.get('search', 'precursor accuracy unit') == 'ppm'

    shifts_and_pime = get_shifts_and_pime(settings)
    
    count = 0
    for result in results:
        mz = get_precursor_mz(result['spectrum'])
        dm_l = acc_l * mz / 1.0e6 if rel else acc_l
        dm_r = acc_r * mz / 1.0e6 if rel else acc_r
        count += 1
        c = result['candidates']
        c = c[c['score'] > score_threshold]
        if min_matched:
            mask = np.array([
                c_[4]['match'] is not None and
                sum(m.sum() for m in c_[4]['match'].values()) >= min_matched
                for c_ in c], dtype=bool)
            c = c[mask]
        mask = []
        for c_ in c:
            mask.append(any(-dm_l < (c_[4]['mzdiff']['Da'] - sh_) / c_[3] < dm_r for sh_ in shifts_and_pime))
        c = c[np.array(mask, dtype=bool)]
    
        if (not c.size) and not show_empty:
            continue
        result['candidates'] = c[:num_candidates]
        yield result
    print 'Unfiltered results:', count
    
def get_shifts_and_pime(settings):
    pime = settings.getint('search', 'precursor isotope mass error') 
    shifts =[float(x) for x in settings.get('search', 'shifts').split(',')]
    dM = mass.nist_mass['C'][13][0] - mass.nist_mass['C'][12][0]
    shifts_and_pime = shifts[:]
    for i in range(pime):
        shifts_and_pime += [x + (i + 1) * dM for x in shifts]
    return shifts_and_pime
    
def write_pepxml(inputfile, settings, results):
    from lxml import etree
    from time import strftime
    from os import path

    if settings.has_option('output', 'path'):
        outpath = settings.get('output', 'path')
    else:
        outpath = path.dirname(inputfile)

    db = settings.get('input', 'database')
    add_decoy = settings.getboolean('input', 'add decoy')
    prefix = settings.get('input', 'decoy prefix')
    mode = settings.get('input', 'decoy method')
    if add_decoy:
        ft = tempfile.NamedTemporaryFile(mode='w')
        fasta.write_decoy_db(db, ft, mode=mode, prefix=prefix)
        settings.set('input', 'database', ft.name)
        settings.set('input', 'add decoy', 'no')

    filename = path.join(outpath, path.splitext(path.basename(inputfile))[0] + path.extsep + 'pep' + path.extsep + 'xml')
    enzyme = settings.get('search', 'enzyme')
    mc = settings.getint('search', 'number of missed cleavages')
    minlen = settings.getint('search', 'peptide minimum length')
    maxlen = settings.getint('search', 'peptide maximum length')
    prefix = settings.get('input', 'decoy prefix')
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
    seen_target.clear()
    seen_decoy.clear()
    for desc, prot in prot_gen(settings):
        dbinfo = desc.split(' ')[0]
        prots[dbinfo] = desc
        for pep in prot_peptides(prot, get_enzyme(enzyme), mc, minlen, maxlen, desc.startswith(prefix)):
            if pep in peptides:
                pept_prot.setdefault(pep, []).append(dbinfo)

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
                vmods.add(aa + k)

    leg = {}
    if settings.has_option('misc', 'legend'):
        leg = settings.get('misc', 'legend')

    ntermcleavage = settings.getfloat('modifications', 'protein nterm cleavage')
    ctermcleavage = settings.getfloat('modifications', 'protein cterm cleavage')

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
                match = candidate[4]['match']
                if match is None: break
                tmp3 = etree.Element('search_hit')
                tmp3.set('hit_rank', str(i + 1))
                mod_sequence = str(candidate[1])
                mod_sequence = normalize_mods(mod_sequence, settings)

                sequence = re.sub(r'[^A-Z]', '', mod_sequence)
                if sequence not in pept_prot:
                    flag = 0
                    print 'WTF'
                    print sequence
                    print mod_sequence
                    sys.stdout.flush()
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
                    tmp3.set('num_matched_ions', str(sum(v.sum() for v in match.values())))
                    tmp3.set('tot_num_ions', '7')  # ???
                    neutral_mass_theor = cmass.fast_mass(sequence, aa_mass=aa_mass)
                    tmp3.set('calc_neutral_pep_mass', str(neutral_mass_theor))
                    tmp3.set('massdiff', str(candidate[4]['mzdiff']['Da']))
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

                    try:
                        aalist = parser.parse(mod_sequence)
                    except:
                        aalist = [a[::-1] for a in parser.parse(mod_sequence[::-1])][::-1]
                    tmp4 = etree.Element('modification_info')
                    ntermmod = 0
                    for idx, aminoacid in enumerate(aalist):
                        if aminoacid in fmods or aminoacid in vmods:
                            if aminoacid.endswith('-') and idx == 0:
                                ntermmod = 1
                                tmp4.set('mod_nterm_mass', str(str(aa_mass.get(aminoacid) + ntermcleavage)))
                            elif aminoacid.startswith('-') and idx == len(aalist) - 1:
                                tmp4.set('mod_cterm_mass', str(aa_mass.get(aminoacid) + ctermcleavage))
                            else:
                                tmp5 = etree.Element('mod_aminoacid_mass')
                                tmp5.set('position', str(idx + 1 - ntermmod))
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

                    tmp4 = etree.Element('search_score')
                    tmp4.set('name', 'fragmentMT')
                    tmp4.set('value', str(candidate[6]))
                    tmp3.append(copy(tmp4))

                    tmp2.append(copy(tmp3))
            if flag:
                tmp.append(copy(tmp2))
                child1.append(copy(tmp))

    s = etree.tostring(root, pretty_print=True)
    output.write(s)

    output.close()
