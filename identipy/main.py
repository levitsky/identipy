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
db = '/home/lev/Downloads/uniprot_sprot.fasta'
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

