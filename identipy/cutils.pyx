cimport cython
from cpython.sequence cimport PySequence_GetSlice
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem
from cpython.float cimport PyFloat_AsDouble
from cpython.tuple cimport PyTuple_GetItem

from pyteomics import cmass
from pyteomics.cmass cimport fast_mass
from pyteomics import electrochem as ec
import numpy as np
cimport numpy as np

np.import_array()

cdef np.dtype dtype_t = np.PyArray_DescrFromType(float)
print(dtype_t)


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
cdef tuple ctheor_spectrum(str peptide, float acc_frag, tuple types,
                           int maxcharge, bint reshape, dict kwargs):
    cdef int pl, charge, i, n, i_type, n_types
    cdef bint nterminal
    cdef str ion_type, maxpart, part
    cdef float maxmass, part_mass
    cdef dict peaks, theoretical_set
    cdef dict aa_mass, ion_comp, mass_data
    cdef set theoretical_set_item, ions_scaled
    cdef np.ndarray[float, ndim=1, mode='c'] marr

    aa_mass = kwargs.get("aa_mass")
    ion_comp = kwargs.get("ion_comp")
    mass_data = kwargs.get("mass_data")

    peaks = {}
    theoretical_set = dict()

    pl = len(peptide) - 1
    n_types = len(types)
    for charge in range(1, maxcharge + 1):
        for i_type in range(n_types):
            ion_type = <str>PyTuple_GetItem(types, i_type)
            nterminal = ion_type[0] in 'abc'
            if nterminal:
                maxpart = <str>PySequence_GetSlice(peptide, 0, -1)
                maxmass = fast_mass(
                    maxpart, ion_type, charge, mass_data, aa_mass,
                    ion_comp)
                marr = np.zeros(pl, dtype=float)
                marr[0] = maxmass
                for i in range(1, pl):
                    part = <str>maxpart[-i]
                    part_mass = PyFloat_AsDouble(<object>PyDict_GetItem(aa_mass, part))
                    marr[i] = marr[i - 1] - part_mass / charge
            else:
                maxpart = <str>PySequence_GetSlice(peptide, 1, pl + 2)
                maxmass = fast_mass(
                    maxpart, ion_type, charge, mass_data, aa_mass,
                    ion_comp)
                marr = np.zeros(pl, dtype=float)
                marr[pl-1] = maxmass
                for i in range(pl-2, -1, -1):
                    part = <str>maxpart[-(i + 2)]
                    part_mass = PyFloat_AsDouble(<object>PyDict_GetItem(aa_mass, part))
                    marr[i] = marr[i + 1] - part_mass / charge

            ions_scaled = set()
            for i in range(marr.shape[0]):
                ions_scaled.add(<int>(marr[i] / acc_frag))

            if ion_type in theoretical_set:
                theoretical_set_item = <set>PyDict_GetItem(theoretical_set, ion_type)
                theoretical_set_item.update(ions_scaled)
            else:
                theoretical_set[ion_type] = ions_scaled

            if not reshape:
                marr.sort()
            else:
                n = marr.size
                marr = marr.reshape((n, 1))
            peaks[ion_type, charge] = marr
    return peaks, theoretical_set


def theor_spectrum(peptide, acc_frag, types=('b', 'y'), maxcharge=None, reshape=False, **kwargs):
    if not maxcharge:
        maxcharge = 1 + int(ec.charge(peptide, pH=2))
    return ctheor_spectrum(peptide, acc_frag, tuple(types), maxcharge, reshape, kwargs)