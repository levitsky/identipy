import numpy as np
from scipy.spatial import cKDTree
from pyteomics import parser, mass, fasta
from itertools import chain
import re
import os
from math import factorial
from tempfile import gettempdir
import sys
from . import scoring, utils 
from multiprocessing import Queue, Process, cpu_count
db = 'uniprot_sprot.fasta'
MAXLEN = 30
acc = 0.02
MC = 1
enzyme = parser.expasy_rules['trypsin']
MINMASS = 300
MAXMASS = 1500

def top_candidates(spectrum, score, seqs, masses, n=1):
    exp_mass = utils.neutral_masses(spectrum)
    candidates = []
    for m, c in exp_mass:
        start = masses.searchsorted(m - acc)
        end = masses.searchsorted(m + acc)
        candidates.extend(seqs[start:end])
    spectrum['__KDTree'] = cKDTree(spectrum['m/z array'].reshape(
        (spectrum['m/z array'].size, 1)))
    return sorted(((score(spectrum, x), x) for x in candidates),
            reverse=True)[:n]

def generate_arrays(folder=gettempdir()):
    if not os.path.isfile(os.path.join(folder, 'seqs.npy')):
        seqs = np.fromiter((pep for _, prot in fasta.read(db)
            for pep in parser.cleave(prot, enzyme, MC)
            if len(pep) < MAXLEN and parser.valid(pep)),
            dtype=np.dtype((np.str_, MAXLEN)))
        masses = np.empty(seqs.shape, dtype=np.float32)
        for i in np.arange(seqs.size):
            masses[i] = mass.fast_mass(seqs[i])
        seqs, ind = np.unique(seqs, return_index=True)
        masses = masses[ind]
        idx = np.argsort(masses)
        masses = masses[idx]
        seqs = seqs[idx]
        np.save(os.path.join(folder, 'masses.npy'), masses)
        np.save(os.path.join(folder, 'seqs.npy'), seqs)
    else:
        seqs = np.load(os.path.join(folder, 'seqs.npy'))
        masses = np.load(os.path.join(folder, 'masses.npy'))
    return masses, seqs

def process_file(f, score, top=1):
    masses, seqs = generate_arrays()
    return ((s, top_candidates(s, score, seqs, masses, top))
            for s in f)

def make_worker(masses, seqs, score, top):
    def worker(qin, qout):
        for spectrum in iter(qin.get, None):
            result = top_candidates(spectrum, score, seqs, masses, top)
            qout.put((spectrum, result))
    return worker

def process_parallel(f, score, top=1):
    masses, seqs = generate_arrays()
    qin = Queue()
    qout = Queue()
    N = cpu_count()
    count = 0
    global worker
    worker = make_worker(masses, seqs, score, top)
    for _ in range(N):
        Process(target=worker, args=(qin, qout)).start()
    for s in f:
        qin.put(s)
        count += 1
    for _ in range(N):
        qin.put(None)
    while count:
        yield qout.get()
        count -= 1

