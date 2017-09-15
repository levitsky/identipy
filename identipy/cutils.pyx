from pyteomics import cmass
import re
from pyteomics import mass, electrochem as ec, auxiliary as aux, fasta, mzml, parser
import numpy as np
import tempfile
import os

def theor_spectrum(peptide, acc_frag, types=('b', 'y'), maxcharge=None, reshape=False, **kwargs):
    cdef int pl, charge, i, n
    # cdef int i, cl
    cdef str ion_type, maxpart
    cdef float maxmass
    peaks = {}
    theoretical_set = dict()#defaultdict(set)#set()
    # theoretical_set = set()
    pl = len(peptide) - 1
    if not maxcharge:
        maxcharge = 1 + int(ec.charge(peptide, pH=2))
    for charge in range(1, maxcharge + 1):
        for ion_type in types:
            nterminal = ion_type[0] in 'abc'
            if nterminal:
                maxpart = peptide[:-1]
                maxmass = cmass.fast_mass(maxpart, ion_type=ion_type, charge=charge, **kwargs)
                marr = np.zeros((pl, ), dtype=float)
                marr[0] = maxmass
                for i in range(1, pl):
                    marr[i] = marr[i-1] - kwargs['aa_mass'][maxpart[-i]]/charge
            else:
                maxpart = peptide[1:]
                maxmass = cmass.fast_mass(maxpart, ion_type=ion_type, charge=charge, **kwargs)
                marr = np.zeros((pl, ), dtype=float)
                marr[pl-1] = maxmass
                for i in range(pl-2, -1, -1):
                    marr[i] = marr[i+1] - kwargs['aa_mass'][maxpart[-(i+2)]]/charge

            tmp = marr / acc_frag
            # tmp = tmp.astype(int)
            # tmp = np.concatenate((tmp, tmp-1, tmp+1))
            # theoretical_set.update(tmp)
            if ion_type in theoretical_set:
                theoretical_set[ion_type].update(set(tmp.astype(int)))
            else:
                theoretical_set[ion_type] = set(tmp.astype(int))

            if not reshape:
                marr.sort()
            else:
                n = marr.size
                marr = marr.reshape((n, 1))
            peaks[ion_type, charge] = marr
    return peaks, theoretical_set