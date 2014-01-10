import re
from pyteomics import mass, electrochem as ec, auxiliary as aux, parser, fasta
import sys
import numpy as np
from multiprocessing import Queue, Process, cpu_count
from string import punctuation
from copy import copy


def find_nearest(array, value):
    return (np.abs(array - value)).argmin()


def get_info(spectrum, result, settings, aa_mass=None):
    'Returns neutral mass, charge state and retention time of the spectrum'
    if not aa_mass:
        aa_mass = get_aa_mass(settings)
    if 'params' in spectrum:
        exp_mass = spectrum['params']['pepmass'][0]
        charge = spectrum['params'].get('charge')
        RT = spectrum['params'].get('rtinseconds')
    else:
        ion = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
        charge = int(ion['charge state'])
        exp_mass = ion['selected ion m/z']
        RT = spectrum['scanList']['scan'][0]['scan start time']
    if isinstance(charge, str):
        states = aux._parse_charge(charge)
        if isinstance(states, int):
            states = [states]
    elif charge is None:
        states = list(range(settings.getint('search', 'minimum charge'), 1 + settings.getint('search', 'maximum charge')))
    else:
        states = [charge]
    states = [s for s in states if settings.getint('search', 'minimum charge') <= s <= settings.getint('search', 'maximum charge')]
    states.sort()
    masses = np.array([exp_mass * ch - ch * mass.nist_mass['H+'][0][0]
            for ch in states])
    idx = find_nearest(masses, mass.fast_mass(result['candidates'][0][1], aa_mass=aa_mass))
    return (masses[idx], states[idx], RT)


def decode(func):
    if sys.version_info.major == 3:
        def f(s, *args, **kwargs):
            if isinstance(s, bytes):
                return func(s.decode('ascii'), *args, **kwargs)
            else:
                return func(s, *args, **kwargs)
        return f
    return func


@decode
def theor_spectrum(peptide, types=('b', 'y'), maxcharge=None, **kwargs):
    peaks = {}
    if not maxcharge:
        maxcharge = 1 + int(ec.charge(peptide, pH=2))
    for ion_type in types:
        ms = []
        for i in range(1, len(peptide) - 1):
            for charge in range(1, maxcharge + 1):
                if ion_type[0] in 'abc':
                    ms.append(mass.fast_mass(
                        str(peptide)[:i], ion_type=ion_type, charge=charge,
                        **kwargs))
                else:
                    ms.append(mass.fast_mass(
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
        ion = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
        charge = int(ion['charge state'])
        exp_mass = ion['selected ion m/z']
    if isinstance(charge, str):
        states = aux._parse_charge(charge)
        if isinstance(states, int):
            states = [states]
    elif charge is None:
        states = list(range(settings.getint('search', 'minimum charge'), 1 + settings.getint('search', 'maximum charge')))
    else:
        states = [charge]
    states = [s for s in states if settings.getint('search', 'minimum charge') <= s <= settings.getint('search', 'maximum charge')]
    states.sort()
    return zip((exp_mass * ch - ch * mass.nist_mass['H+'][0][0]
            for ch in states), states)


@aux.memoize(10)
def import_(name):
    """Import a function by name: module.function or
    module.submodule.function, etc. Return the function object."""

    mod, f = name.rsplit('.', 1)
    return getattr(__import__(mod, fromlist=[f]), f)


def get_aa_mass(settings):
    if settings.has_option('misc', 'aa_mass'):
        return settings.get('misc', 'aa_mass')
    aa_mass = mass.std_aa_mass.copy()
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

    return aa_mass


def multimap(n, func, it):
    if n == 0:
        try:
            n = cpu_count()
        except NotImplementedError:
            n = 1
    if n == 1:
        for s in it:
            yield func(s)
    else:
        def worker(qin, qout):
            for spectrum in iter(qin.get, None):
                result = func(spectrum)
                qout.put(result)
        qin = Queue()
        qout = Queue()
        count = 0
        for _ in range(n):
            Process(target=worker, args=(qin, qout)).start()
        for s in it:
            qin.put(s)
            count += 1
        for _ in range(n):
            qin.put(None)
        while count:
            yield qout.get()
            count -= 1


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


def write_pepxml(inputfile, settings, results):
    from lxml import etree
    from time import strftime
    from os import path

    if settings.has_option('misc', 'aa_mass'):
        aa_mass = settings.get('misc', 'aa_mass')
    else:
        aa_mass = get_aa_mass(settings)

    filename = path.splitext(inputfile)[0] + path.extsep + 'pep' + path.extsep + 'xml'
    enzyme = settings.get('search', 'enzyme')
    search_engine = 'IdentiPy'
    database = settings.get('input', 'database')
    missed_cleavages = settings.getint('search', 'missed cleavages')


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

    results = [x for x in results if x['candidates']]
    pept_prot = dict()
    prots = dict()
    peptides = set(x['candidates'][i][1] for x in results for i in range(
                settings.getint('output', 'candidates') or len(x['candidates'])))
    for desc, prot in fasta.read(database):
        prots[desc.split(' ')[0]] = desc
        for pep in parser.cleave(prot, parser.expasy_rules.get(enzyme, enzyme), missed_cleavages):
            if pep in peptides:
                dbinfo = desc
                try:
                    pept_prot[pep] = np.append(pept_prot[pep], dbinfo.split(' ')[0])
                except:
                    pept_prot[pep] = np.array([dbinfo.split(' ')[0]], dtype = '|S50')

    for idx, result in enumerate(results):
        if result['candidates']:
            tmp = etree.Element('spectrum_query')
            spectrum = result['spectrum']
            tmp.set('spectrum', spectrum['params']['title'])
            tmp.set('start_scan', str(idx))  # ???
            tmp.set('end_scan', str(idx))  # ???
            tmp.set('index', str(idx))  # ???

            neutral_mass, charge_state, RT = get_info(spectrum, result, settings, aa_mass)
            tmp.set('precursor_neutral_mass', str(neutral_mass))
            tmp.set('assumed_charge', str(charge_state))
            if RT:
                tmp.set('retention_time_sec', str(RT))

            tmp2 = etree.Element('search_result')
            result['candidates'] = result['candidates'][:len(result['e-values'])]

            for i, candidate in enumerate(result['candidates']):
                tmp3 = etree.Element('search_hit')
                tmp3.set('hit_rank', str(i + 1))
                sequence = candidate[1]
                tmp3.set('peptide', sequence)
                tmp3.set('peptide_prev_aa', 'K')  # ???
                tmp3.set('peptide_next_aa', 'K')  # ???
                proteins = pept_prot[sequence]

                tmp3.set('protein', prots[proteins[0]].split(' ', 1)[0])
                tmp3.set('protein_descr', prots[proteins[0]].split(' ', 1)[1])

                num_tot_proteins = len(proteins)
                tmp3.set('num_tot_proteins', str(num_tot_proteins))
                tmp3.set('num_matched_ions', '7')  # ???
                tmp3.set('tot_num_ions', '7')  # ???
                neutral_mass_theor = mass.fast_mass(sequence, aa_mass=aa_mass)
                tmp3.set('calc_neutral_pep_mass', str(neutral_mass_theor))
                tmp3.set('massdiff', str(neutral_mass - neutral_mass_theor))
                tmp3.set('num_tol_term', '2')  # ???
                tmp3.set('num_missed_cleavages', str(len(parser.cleave(sequence, parser.expasy_rules[enzyme], 0)) - 1))
                tmp3.set('is_rejected', '0')  # ???

                if num_tot_proteins > 1:
                    for idx in range(len(proteins)):
                        if idx != 0:
                            tmp4 = etree.Element('alternative_protein')
                            tmp4.set('protein', prots[proteins[idx]].split(' ', 1)[0])
                            tmp4.set('protein_descr', prots[proteins[idx]].split(' ', 1)[1])
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

                tmp2.append(copy(tmp3))
            tmp.append(copy(tmp2))
            child1.append(copy(tmp))

    s = etree.tostring(root, pretty_print=True)
    output.write(s)

    output.close()
