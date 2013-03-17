from __future__ import print_function
import numpy as np
from scipy.spatial import cKDTree
from pyteomics import parser, mass, fasta, auxiliary as aux, electrochem as ec, mgf
from itertools import chain
import re
import os
from math import factorial
from tempfile import gettempdir
import sys
db = '/home/lev/Downloads/uniprot_sprot.fasta'
MAXLEN = 30
acc = 0.02
MC = 1
enzyme = parser.expasy_rules['trypsin']
MINMASS = 300
MAXMASS = 1500

def decode(func):
    if sys.version_info.major == 3:
        def f(s, *args, **kwargs):
            if isinstance(s, bytes):
                return func(s.decode('ascii'), *args, **kwargs)
        return f
    return func

@decode
def theor_spectrum(peptide, types=('y', 'b'), maxcharge=None):
    peaks = {}
    if not maxcharge:
        maxcharge = 1 + int(ec.charge(peptide, pH=2))
    for ion_type in types:
        ms = []
        for i in range(1, len(peptide)-1):
            for charge in range(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    ms.append(mass.fast_mass(
                        str(peptide)[:i], ion_type=ion_type, charge=charge))
                else:
                    ms.append(mass.fast_mass(
                        str(peptide)[i:], ion_type=ion_type, charge=charge))
        marr = np.array(ms)
        marr.sort()
        peaks[ion_type] = marr
    return peaks

def simple_score(spectrum, peptide):
    charge = max(c for m, c in neutral_masses(spectrum))
    theor = theor_spectrum(peptide, maxcharge=charge)
    fragments = np.concatenate(theor.values())
    dist, ind = spectrum['__KDTree'].query(fragments.reshape((fragments.size, 1)),
            distance_upper_bound = acc)
    return spectrum['intensity array'][ind[dist != np.inf]].sum()

def hyperscore(spectrum, peptide):
    """A simple implementation of X!Tandem's Hyperscore."""

    charge = max(c for m, c in neutral_masses(spectrum))
    theor = theor_spectrum(peptide, maxcharge=charge)
    score = 0
    mult = []
    for fragments in theor.values():
        n = fragments.size
        dist, ind = spectrum['__KDTree'].query(fragments.reshape((n, 1)),
            distance_upper_bound=acc, eps=acc)
        mask = (dist != np.inf)
        mult.append(mask.sum())
        score += spectrum['intensity array'][ind[mask]].sum()
    if score:
        for m in mult: score *= m
    return score

def top_candidates(spectrum, score, seqs, masses, n=1):
    exp_mass = neutral_masses(spectrum)
    candidates = []
    for m, c in exp_mass:
        start = masses.searchsorted(m - acc)
        end = masses.searchsorted(m + acc)
        candidates.extend(seqs[start:end])
    spectrum['__KDTree'] = cKDTree(spectrum['m/z array'].reshape(
        (spectrum['m/z array'].size, 1)))
    return sorted(((score(spectrum, x), x) for x in candidates),
            reverse=True)[:n]

def neutral_masses(spectrum):
    if 'params' in spectrum:
        exp_mass = spectrum['params']['pepmass'][0]
        charge = spectrum['params'].get('charge')
    else:
        exp_mass = spectrum['base peak']
        charge = spectrum['charge']
    if isinstance(charge, str):
        states = re.findall(r'(\d+)(\+|-)?', charge)
        states = [int(num)*(-1 if sign == '-' else 1)
                for num, sign in states]
    else: states = [charge]
    states.sort()
    return zip((exp_mass*ch - ch*mass.nist_mass['H+'][0][0]
            for ch in states), states)

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

if __name__ == '__main__':
    with mgf.read(os.path.join(os.pardir, 'tests', 'swedcad.mgf')) as r:
        true = 0
        for spectrum, result in process_file(r, hyperscore, 1):
            if result and (result[0][1] == spectrum['params']['title']):
                true += 1
            else:
                print('{} != {}'.format(result[0][1], spectrum['params']['title']))
    print('Correct:', true)
