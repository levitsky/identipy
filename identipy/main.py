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
from ConfigParser import ConfigParser
from . import scoring, utils 

def top_candidates_from_arrays(spectrum, settings,
        score=None, seqs=None, masses=None, acc=None, rel=False, n=1):
    # todo: handle unspecified parameters (exctract from settings)
    exp_mass = utils.neutral_masses(spectrum)
    candidates = []
    for m, c in exp_mass:
        dm = acc*c*m/1.0e6 if rel else acc*c
        start = masses.searchsorted(m - dm)
        end = masses.searchsorted(m + dm)
        candidates.extend(seqs[start:end])
    spectrum['__KDTree'] = cKDTree(spectrum['m/z array'].reshape(
        (spectrum['m/z array'].size, 1)))
    return sorted(((score(spectrum, x, settings), x) for x in candidates),
            reverse=True)[:n]

def get_arrays(settings):
    db = settings.get('input', 'database')
    folder = settings.get('performance', 'folder')
    enzyme = settings.get('search', 'enzyme')
    enzyme = parser.expasy_rules.get(enzyme, enzyme)
    mc = settings.getint('search', 'miscleavages')
    minlen = settings.getint('search', 'peptide minimum length')
    maxlen = settings.getint('search', 'peptide maximum length')
    index = os.path.join(folder, 'identipy.idx')

    profile = (db, enzyme, mc, minlen, maxlen)
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
            if minlen <= len(pep) <= maxlen and parser.valid(pep)),
            dtype=np.dtype((np.str_, maxlen)))
        masses = np.empty(seqs.shape, dtype=np.float32)
        for i in np.arange(seqs.size):
            masses[i] = mass.fast_mass(seqs[i])
        seqs, ind = np.unique(seqs, return_index=True)
        masses = masses[ind]
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

def process_file(f, settings):
    # prepare the function
    mode = settings.get('performance', 'pre-calculation')
    score = getattr(scoring, settings.get('scoring', 'score'))
    prec_acc = settings.getfloat('search', 'precursor accuracy value')
    unit = settings.get('search', 'precursor accuracy unit')
    if unit == 'ppm':
        rel = True
    elif unit in {'Th', 'Da', 'amu'}:
        rel = False
    else:
        raise ValueError('Unrecognized precursor accuracy unit')
    
    if mode == 'some': # work with numpy arrays
        masses, seqs = get_arrays(settings)
        func = lambda s: top_candidates_from_arrays(s, settings, score,
                seqs, masses,
                prec_acc, rel, settings.getint('output', 'candidates'))
    else:
        raise NotImplementedError

    # decide on multiprocessing
    n = settings.getint('performance', 'processes')
    if n == 0:
        try:
            n = cpu_count()
        except NotImplementedError:
            n = 1
    if n == 1:
        for s in f:
            yield s, func(s)
    else:
        def worker(qin, qout):
            for spectrum in iter(qin.get, None):
                result = func(spectrum)
                qout.put((spectrum, result))
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
    config = ConfigParser()
    if default_name:
        config.read(default_name)
    if fname:
        config.read(fname)
    return config
        

#def make_worker(masses, seqs, score, top):
#    def worker(qin, qout):
#        for spectrum in iter(qin.get, None):
#            result = top_candidates(spectrum, score, seqs, masses, top)
#            qout.put((spectrum, result))
#    return worker
#
#def process_parallel(f, settings):
#    masses, seqs = generate_arrays()
#    qin = Queue()
#    qout = Queue()
#    N = cpu_count()
#    count = 0
#    global worker
#    worker = make_worker(masses, seqs, score, top)
#    for _ in range(N):
#        Process(target=worker, args=(qin, qout)).start()
#    for s in f:
#        qin.put(s)
#        count += 1
#    for _ in range(N):
#        qin.put(None)
#    while count:
#        yield qout.get()
#        count -= 1

