import numpy as np
from pyteomics import parser, mass, fasta, auxiliary as aux, mgf, mzml
from itertools import chain
import re
import os
from math import factorial
from tempfile import gettempdir, NamedTemporaryFile
import ast
import sys
import hashlib
from copy import copy
from string import punctuation
from . import scoring, utils 
try:
    from configparser import RawConfigParser
except ImportError:
    from ConfigParser import RawConfigParser

def top_candidates_from_arrays(spectrum, settings):
    dynrange = settings.getfloat('scoring', 'dynamic range')
    if dynrange:
        i = spectrum['intensity array'] > spectrum['intensity array'].max(
                ) / dynrange
        spectrum['intensity array'] = spectrum['intensity array'][i]
        spectrum['m/z array'] = spectrum['m/z array'][i]
    maxpeaks = settings.getint('scoring', 'maximum peaks')
    minpeaks = settings.getint('scoring', 'minimum peaks')
    if minpeaks > maxpeaks:
        raise ValueError('minpeaks > maxpeaks: {} and {}'.format(
            minpeaks, maxpeaks))
    if minpeaks and spectrum['intensity array'].size < minpeaks:
        return []
    if maxpeaks and spectrum['intensity array'].size > maxpeaks:
        i = np.argsort(spectrum['intensity array'])[-maxpeaks:]
        spectrum['intensity array'] = spectrum['intensity array'][i]
        spectrum['m/z array'] = spectrum['m/z array'][i]

    masses, seqs = get_arrays(settings) if not settings.has_option(
            'performance', 'arrays') else settings.get('performance', 'arrays')
    exp_mass = utils.neutral_masses(spectrum)
    n = settings.getint('output', 'candidates')
    score = utils.import_(settings.get('scoring', 'score'))
    acc = settings.getfloat('search', 'precursor accuracy value')
    unit = settings.get('search', 'precursor accuracy unit')
    if unit == 'ppm':
        rel = True
    elif unit in {'Th', 'Da', 'amu'}:
        rel = False
    else:
        raise ValueError('Unrecognized precursor accuracy unit: ' + unit)
 
    candidates = []
    for m, c in exp_mass:
        dm = acc*m/1.0e6 if rel else acc*c
        start = masses.searchsorted(m - dm)
        end = masses.searchsorted(m + dm)
        candidates.extend(seqs[start:end])

    threshold = settings.getfloat('scoring', 'score threshold')
    result = [(score(spectrum, x, settings), x) for x in candidates]
    result = sorted((x for x in result if x[0] > threshold), reverse=True)[:n]
    if settings.has_option('misc', 'legend'):
        mods = list(zip(settings.get('misc', 'legend'), punctuation))
        res = []
        for score, cand in result:
            for (mod, aa), char in mods:
                cand = cand.replace(char, mod+aa)
            res.append((score, cand))
        result = res
    return result

#@aux.memoize(10)
def get_arrays(settings):
    db = settings.get('input', 'database')
    hasher = settings.get('misc', 'hash')
    dbhash = hashlib.new(hasher)
    with open(db) as f:
        for line in f:
            dbhash.update(line.encode('ascii'))
    dbhash = dbhash.hexdigest()
    folder = settings.get('performance', 'folder')
    enzyme = settings.get('search', 'enzyme')
    enzyme = parser.expasy_rules.get(enzyme, enzyme)
    mc = settings.getint('search', 'miscleavages')
    minlen = settings.getint('search', 'peptide minimum length')
    maxlen = settings.getint('search', 'peptide maximum length')
    aa_mass = utils.aa_mass(settings)
    add_decoy = settings.getboolean('input', 'add decoy')
    index = os.path.join(folder, 'identipy.idx')

    profile = (dbhash, add_decoy, enzyme, mc, minlen, maxlen, aa_mass)
    arr_name = None
    if os.path.isfile(index):
        with open(index) as ifile:
            for line in ifile:
                t, name = line.strip().split('\t')
                if ast.literal_eval(t) == profile:
                    arr_name = name
                    break

    if arr_name is None:
        def peps():
            if not add_decoy:
                prots = (prot for _, prot in fasta.read(db))
            else:
                prefix = settings.get('input', 'decoy prefix')
                mode = settings.get('input', 'decoy method')
                prots = (prot for _, prot in fasta.decoy_db(db, mode=mode,
                    prefix=prefix))
            func = lambda prot: [pep for pep in parser.cleave(prot, enzyme, mc)
                    if minlen <= len(pep) <= maxlen and parser.fast_valid(pep)]
            n = settings.getint('performance', 'processes')
            return chain.from_iterable(
                    utils.multimap(n, func, prots))

        seqs = np.fromiter(peps(), dtype=np.dtype((np.str_, maxlen)))
        seqs = np.unique(seqs)
        masses = np.empty(seqs.shape, dtype=np.float32)
        for i in np.arange(seqs.size):
            masses[i] = mass.fast_mass(seqs[i], aa_mass=aa_mass)
        idx = np.argsort(masses)
        masses = masses[idx]
        seqs = seqs[idx]
        with NamedTemporaryFile(suffix='.npz', prefix='identipy_', dir=folder,
                delete=False) as outfile:
            np.savez_compressed(outfile, masses=masses, seqs=seqs)
            name = outfile.name
        with open(index, 'a') as ifile:
            ifile.write(str(profile) + '\t' + name + '\n')

    else:
        npz = np.load(arr_name)
        seqs = npz['seqs']
        masses = npz['masses']
    return masses, seqs

def spectrum_processor(settings):
    processor = settings.get('misc', 'spectrum processor')
    if '.' in processor:
        return utils.import_(processor)(settings)
   
    mode = settings.get('performance', 'pre-calculation')
    if mode == 'some': # work with numpy arrays
        settings = copy(settings)
        settings.set('performance', 'arrays', get_arrays(settings))
        candidates = lambda s: top_candidates_from_arrays(s, settings)
        if processor == 'minimal':
            return lambda s: {'spectrum': s, 'candidates': candidates(s)}
        elif processor == 'e-value':
            def f(s):
                c = candidates(s)
                return  {'spectrum': s, 'candidates': c,
                    'e-values': scoring.evalues(c)}
            return f
    else:
        raise NotImplementedError

def process_spectra(f, settings):
    # prepare the function
    func = spectrum_processor(settings)

    # decide on multiprocessing
    n = settings.getint('performance', 'processes')
    return utils.multimap(n, func, f)

def process_file(fname, settings):
    stage1 = settings.get('misc', 'first stage')
    if stage1:
        return double_run(fname, settings, utils.import_(stage1))
    if settings.get('modifications', 'variable'):
        return double_run(fname, settings, varmod_stage1)
    else:
        return process_spectra({'mgf': mgf, 'mzml': mzml
            }[fname.rsplit('.', 1)[-1].lower()].read(fname), settings)

def double_run(fname, settings, stage1):
    new_settings = stage1(fname, settings)
    return process_file(fname, new_settings)

def varmod_stage1(fname, settings):
    """Take mods, make a function that yields new settings"""
    mods = settings.get('modifications', 'variable')
    mods = [parser._split_label(l) for l in re.split(r',\s*', mods)]
    mods.sort(key=lambda x: len(x[0]), reverse=True)
    assert all(len(m) == 2 for m in mods), 'unmodified residue given'
    mod_dict = {}
    for mod, aa in mods:
        mod_dict.setdefault(mod, []).append(aa)

    candidates = set()
    settings = copy(settings)
    settings.set('modifications', 'variable', '')
    for s in process_file(fname, settings):
        candidates.update(c[1] for c in s['candidates'])

    n = settings.getint('modifications', 'maximum variable mods')
    seq_iter = chain.from_iterable(
            (parser.tostring(x, False) for x in 
                parser.isoforms(seq, variable_mods=mod_dict, format='split')
                if ((len(x[0]) > 2) + (len(x[-1]) > 2) + sum(
                    len(y) > 1 for y in x[1:-1])) <= n)
            for seq in candidates)
    def prepare_seqs():
        for seq in seq_iter:
            for (mod, aa), char in zip(mods, punctuation):
                seq = seq.replace(mod+aa, char)
            yield seq
    maxlen = settings.getint('search', 'peptide maximum length')
    seqs = np.fromiter(prepare_seqs(), dtype=np.dtype((np.str_, maxlen)))
    aa_mass = mass.std_aa_mass.copy()
    for (mod, aa), char in zip(mods, punctuation):
        aa_mass[char] = aa_mass[aa] + settings.getfloat('modifications', mod)
    masses = np.empty(seqs.shape, dtype=np.float32)
    for i in range(seqs.size):
        masses[i] = mass.fast_mass(seqs[i], aa_mass=aa_mass)
    i = masses.argsort()
    masses = masses[i]
    seqs = seqs[i]
    new_settings = copy(settings)
    new_settings.set('misc', 'aa_mass', aa_mass)
    new_settings.set('misc', 'legend', mods)
    new_settings.set('performance', 'arrays', (masses, seqs))
    maxlen = settings.getint('search', 'peptide maximum length')
    return new_settings

@aux.memoize(10)
def settings(fname=None, default_name=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), os.pardir, 'default.cfg')):
    """Read a configuration file and return a :py:class:`utils.Config` object.
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
#   config = utils.Config()
#   conv = {True: utils.import_, False: ast.literal_eval}
#   for section in raw_config.sections():
#       for option in raw_config.options(section):
#           val = raw_config.get(section, option)
#           config[section][option] = conv[
#                   raw_config.has_option('import', section + ', ' + option)](
#                           raw_config.get(section, option)) if val else None
#   return config
