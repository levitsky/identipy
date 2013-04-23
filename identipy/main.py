import numpy as np
from pyteomics import parser, mass, fasta, auxiliary as aux
from itertools import chain
import re
import os
from math import factorial
from tempfile import gettempdir, NamedTemporaryFile
import ast
import sys
import hashlib
from . import scoring, utils 

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

    masses, seqs = get_arrays(settings)
    exp_mass = utils.neutral_masses(spectrum)
    n = settings.getint('output', 'candidates')
    score = _get_score(settings)
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
    result = [(score(spectrum, x, settings), x.decode('ascii')) for x in candidates]
    return sorted((x for x in result if x[0] > threshold), reverse=True)[:n]

@aux.memoize(10)
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

#@aux.memoize(10)
def spectrum_processor(settings):
    processor = settings.get('misc', 'spectrum processor')
    if '.' in processor:
        return utils.import_function(processor)(settings)
   
    mode = settings.get('performance', 'pre-calculation')
    if mode == 'some': # work with numpy arrays
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

def process_file(f, settings):
    # prepare the function
    func = spectrum_processor(settings)

    # decide on multiprocessing
    n = settings.getint('performance', 'processes')
    return utils.multimap(n, func, f)

#@aux.memoize(10)
def settings(fname=None, default_name=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), os.pardir, 'default.cfg')):
    kw = {'inline_comment_prefixes': ('#', ';')
            } if sys.version_info.major == 3 else {}
    kw['dict_type'] = dict
    config = utils.Config(**kw)
    if default_name:
        config.read(default_name)
    if fname:
        config.read(fname)
    return config

#@aux.memoize(10)
def _get_score(settings):
    score_name = settings.get('scoring', 'score')
    if '.' in score_name:
        score = utils.import_function(score_name)
    else:
        score = getattr(scoring, score_name)
    return score
