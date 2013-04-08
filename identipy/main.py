import numpy as np
from scipy.spatial import cKDTree
from pyteomics import parser, mass, fasta
from itertools import chain
import re
import os
from math import factorial
from tempfile import gettempdir, NamedTemporaryFile
import ast
import sys
from multiprocessing import Queue, Process, cpu_count
import hashlib
try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import ConfigParser
from . import scoring, utils 

def top_candidates_from_arrays(spectrum, settings,
        score=None, seqs=None, masses=None, acc=None, rel=False, n=1):
    # todo: handle unspecified parameters (exctract from settings)
    exp_mass = utils.neutral_masses(spectrum)
    candidates = []
    for m, c in exp_mass:
        dm = acc*m/1.0e6 if rel else acc*c
        start = masses.searchsorted(m - dm)
        end = masses.searchsorted(m + dm)
        candidates.extend(seqs[start:end])
    i = np.nonzero(spectrum['intensity array'])
    spectrum['intensity array'] = spectrum['intensity array'][i]
    spectrum['m/z array'] = spectrum['m/z array'][i]
    spectrum['__KDTree'] = cKDTree(spectrum['m/z array'].reshape(
        (spectrum['m/z array'].size, 1)))
    return sorted(((score(spectrum, x, settings), x.decode('ascii')) for x in candidates),
            reverse=True)[:n]

def get_arrays(settings):
    db = settings.get('input', 'database')
    hasher = settings.get('misc', 'hash')
    dbhash = hashlib.new(hasher)
    with open(db) as f:
        for line in f:
            dbhash.update(line)
    dbhash = dbhash.hexdigest()
    folder = settings.get('performance', 'folder')
    enzyme = settings.get('search', 'enzyme')
    enzyme = parser.expasy_rules.get(enzyme, enzyme)
    mc = settings.getint('search', 'miscleavages')
    minlen = settings.getint('search', 'peptide minimum length')
    maxlen = settings.getint('search', 'peptide maximum length')
    fmods = settings.get('modifications', 'fixed')
    aa_mass = utils._aa_mass(fmods)
    index = os.path.join(folder, 'identipy.idx')

    profile = (dbhash, enzyme, mc, minlen, maxlen, fmods)
    arr_name = None
    if os.path.isfile(index):
        with open(index) as ifile:
            for line in ifile:
                t, name = line.strip().split('\t')
                if ast.literal_eval(t) == profile:
                    arr_name = name
                    break

    if arr_name is None:
        seqs = np.fromiter((pep for _, prot in fasta.read(db)
            for pep in parser.cleave(prot, enzyme, mc)
            if minlen <= len(pep) <= maxlen and parser.fast_valid(pep)),
            dtype=np.dtype((np.str_, maxlen)))
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
        return utils.import_function(processor)(settings)
    mode = settings.get('performance', 'pre-calculation')
    score_name = settings.get('scoring', 'score')
    if '.' in score_name:
        score = utils.import_function(score_name)
    else:
        score = getattr(scoring, score_name)
    prec_acc = settings.getfloat('search', 'precursor accuracy value')
    unit = settings.get('search', 'precursor accuracy unit')
    if unit == 'ppm':
        rel = True
    elif unit in {'Th', 'Da', 'amu'}:
        rel = False
    else:
        raise ValueError('Unrecognized precursor accuracy unit: ' + unit)
    
    if mode == 'some': # work with numpy arrays
        masses, seqs = get_arrays(settings)
        candidates = lambda s: top_candidates_from_arrays(s, settings, score,
                seqs, masses,
                prec_acc, rel, settings.getint('output', 'candidates'))
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
    if n == 0:
        try:
            n = cpu_count()
        except NotImplementedError:
            n = 1
    if n == 1:
        for s in f:
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
        for s in f:
            qin.put(s)
            count += 1
        for _ in range(n):
            qin.put(None)
        while count:
            yield qout.get()
            count -= 1

def settings(fname=None, default_name=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), os.pardir, 'default.cfg')):
    kw = {'inline_comment_prefixes': ('#', ';')
            } if sys.version_info.major == 3 else {}
    kw['dict_type'] = dict
    config = ConfigParser(**kw)
    if default_name:
        config.read(default_name)
    if fname:
        config.read(fname)
    return config
