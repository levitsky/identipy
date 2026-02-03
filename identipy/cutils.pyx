cimport cython
from cpython.sequence cimport PySequence_GetSlice
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem
from cpython.float cimport PyFloat_AsDouble
from cpython.tuple cimport PyTuple_GetItem

from pyteomics import cmass
from math import factorial

cimport pyteomics.cmass as cmass

from pyteomics import electrochem as ec
import numpy as np
cimport numpy as np

np.import_array()


cdef:
    dict std_aa_mass = cmass.std_aa_mass
    dict std_ion_comp = cmass.std_ion_comp
    dict nist_mass = cmass._nist_mass

cdef dict ion_shift_dict

ion_shift_dict = {
    'a': 46.00547930326002,
    'b': 18.010564683699954,
    'c': 0.984015582689949,
    'x': -25.979264555419945,
    'y': 0.0,
    'z': 17.026549101010005,
}


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
def cdeisotope(np.ndarray mz, np.ndarray intens, float acc, int charge, int charge_max=False):

    cdef int i, j
    cdef set skip
    cdef list add, c_range, ix
    cdef float h, c13, search_limit, d

    h = 1.0072765
    c13 = 1.00335
    i = mz.size-2
    skip = set()
    add = []
    c_range = list(range(1, charge+1))
    search_limit = c13 + acc * 1.1

    while i >= 0:
        j = min(mz.size-1, mz.searchsorted(mz[i] + search_limit, side='right'))
        while j > i:
            if intens[i] > intens[j]:
                d = mz[j] - mz[i]
                if d > search_limit:
                    j -= 1
                    continue
                for z in c_range:
                    if abs(d - c13/z) < acc:
                        skip.add(j)
                        if z > 1:
                            add.append((i, z))
                        break
            j -= 1
        i -= 1
    for i, z in add:
        mz[i] = mz[i]*z - (z-1)*h
    if len(skip):
        ix = list(set(range(mz.size)).difference(skip))
        mz = mz[ix]
        intens = intens[ix]
    return mz, intens



@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
def RNHS_ultrafast(dict cur_idict, dict theoretical_set, int min_matched, dict best_res, set allowed_idx, int max_v):

    cdef int total_matched, nm_key, num_b_ions, num_y_ions
    cdef dict cnt_b, cnt_y
    cdef set out
    cdef float best_res_val
 
    if not cur_idict:
        return None

    total_matched = 0

    cnt_b = dict()
    cnt_y = dict()

    for ion in theoretical_set['b']:
        if ion in cur_idict:
            for xx in cur_idict[ion]:
                if xx not in cnt_b:
                    cnt_b[xx] = 1
                else:
                    cnt_b[xx] += 1
            total_matched += 1

    for ion in theoretical_set['y']:
        if ion in cur_idict:
            for xx in cur_idict[ion]:
                if xx not in cnt_y:
                    cnt_y[xx] = 1
                else:
                    cnt_y[xx] += 1
            total_matched += 1

    if total_matched < min_matched:
        return None

    out = set()
    for k in allowed_idx:
        num_b_ions = 0
        num_y_ions = 0
        if k in cnt_b:
            num_b_ions = cnt_b[k]
        if k in cnt_y:
            num_y_ions = cnt_y[k]
        if num_b_ions + num_y_ions >= min_matched:
            best_res_val = best_res.get(k, 0)
            if not best_res_val or -factorial(num_b_ions) * factorial(num_y_ions) <= best_res_val:
                out.add(k)
    return out

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
def RNHS_fast_old(set spectrum_fastset, dict spectrum_idict , dict theoretical_set, int min_matched):
    cdef int matched_approx_b, matched_approx_y, matched_approx
    cdef set matched_b, matched_y
    cdef float isum
    isum = 0
    matched_b = spectrum_fastset.intersection(theoretical_set['b'])
    matched_y = spectrum_fastset.intersection(theoretical_set['y'])
    matched_approx_b = len(matched_b)
    matched_approx_y = len(matched_y)
    matched_approx = matched_approx_b + matched_approx_y
    if matched_approx >= min_matched:
        for fr in matched_b:
            isum += spectrum_idict[fr]
        for fr in matched_y:
            isum += spectrum_idict[fr]
        return matched_approx, factorial(matched_approx_b) * factorial(matched_approx_y) * isum
    else:
        return 0, 0

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
def RNHS_fast(set spectrum_fastset, dict spectrum_idict , dict theoretical_set, int min_matched, dict rank_map):
    cdef int matched_approx_b, matched_approx_y, matched_approx
    cdef set matched_b, matched_y
    cdef float isum
    cdef list all_matched
    isum = 0

    all_matched = []
    for ion in theoretical_set:
        matched_tmp = spectrum_fastset.intersection(theoretical_set[ion])
        all_matched.append((ion, matched_tmp))
    matched_approx = sum(len(z) for z in all_matched)
    if matched_approx >= min_matched:
        for ion, matched_tmp in all_matched:
            for fr in matched_tmp:
                i_rank = spectrum_idict[fr]
                if i_rank in rank_map:
                    tmp_d = rank_map[i_rank]
                    if ion in tmp_d:
                        isum += tmp_d[ion]
                    else:
                        isum += rank_map['m']
        return matched_approx, isum
    else:
        return 0, 0


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
def RNHS_fast_basic(set spectrum_fastset, dict spectrum_idict , dict theoretical_set, int min_matched):
    cdef int matched_approx_b, matched_approx_y, matched_approx
    cdef set matched_b, matched_y
    cdef float isum
    cdef list all_matched
    isum = 0

    all_matched = []
    for ion in theoretical_set:
        matched_tmp = spectrum_fastset.intersection(theoretical_set[ion])
        all_matched.append(matched_tmp)
    matched_approx = sum(len(z) for z in all_matched)
    if matched_approx >= min_matched:
        for matched_tmp in all_matched:
            for fr in matched_tmp:
                isum += spectrum_idict[fr]
        for matched_tmp in all_matched:
            isum *= factorial(len(matched_tmp))
        return matched_approx, isum
    else:
        return 0, 0

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
cdef float calc_ions_from_neutral_mass(str peptide, float nm, str ion_type, int charge, dict aa_mass, float cterm_mass, float nterm_mass):
    cdef float nmi
    if ion_type in 'abc':
        nmi = nm - aa_mass[peptide[-1]] - ion_shift_dict[ion_type] - (cterm_mass - 17.002735)
    else:
        nmi = nm - aa_mass[peptide[0]] - ion_shift_dict[ion_type] - (nterm_mass - 1.007825)
    return (nmi + 1.0072764667700085 * charge) / charge 

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
cdef list get_n_ions(str peptide, float maxmass, int pl, int charge, dict k_aa_mass):
    cdef int i
    cdef list tmp
    tmp = [maxmass, ]
    for i in xrange(1, pl):
        tmp.append(tmp[-1] - k_aa_mass[peptide[-i-1]]/charge)
    return tmp

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
cdef list get_c_ions(str peptide, float maxmass, int pl, int charge, dict k_aa_mass):
    cdef int i
    cdef list tmp
    tmp = [maxmass, ]
    for i in xrange(pl-2, -1, -1):
        tmp.append(tmp[-1] - k_aa_mass[peptide[-(i+2)]]/charge)
    return tmp

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(True)
cdef tuple ctheor_spectrum(str peptide, double acc_frag, double nterm_mass, double cterm_mass, tuple types,
                           int maxcharge, bint reshape, dict kwargs):
    cdef int pl, charge, i, n, i_type, n_types
    cdef bint nterminal
    cdef str ion_type, maxpart, part
    cdef float maxmass, part_mass, nm
    cdef dict peaks, theoretical_set
    cdef dict aa_mass, ion_comp, mass_data
    cdef list theoretical_set_item
    cdef list ions_scaled, marr
    cdef object marr_storage

    aa_mass = kwargs.get("aa_mass")
    if aa_mass is None:
        aa_mass = std_aa_mass
    ion_comp = kwargs.get("ion_comp")
    if ion_comp is None:
        ion_comp = std_ion_comp
    mass_data = kwargs.get("mass_data")
    if mass_data is None:
        mass_data = nist_mass
    nm = kwargs.get("nm")
    if nm is None:
        nm = cmass.fast_mass(peptide, **kwargs) + (nterm_mass - 1.007825) + (cterm_mass - 17.002735)
    peaks = {}
    theoretical_set = dict()

    pl = len(peptide) - 1
    n_types = len(types)
    for charge in range(1, maxcharge + 1):
        for i_type in range(n_types):
            ion_type = <str>PyTuple_GetItem(types, i_type)
            nterminal = ion_type[0] in 'abc'
            if nterminal:
                maxmass = calc_ions_from_neutral_mass(peptide, nm, ion_type=ion_type, charge=charge,
                                aa_mass=kwargs['aa_mass'], cterm_mass=cterm_mass, nterm_mass=nterm_mass)
                marr = get_n_ions(peptide, maxmass, pl, charge, kwargs['aa_mass'])
            else:
                maxmass = calc_ions_from_neutral_mass(peptide, nm, ion_type=ion_type, charge=charge,
                                aa_mass=kwargs['aa_mass'], cterm_mass=cterm_mass, nterm_mass=nterm_mass)
                marr = get_c_ions(peptide, maxmass, pl, charge, kwargs['aa_mass'])

            iname = (ion_type, charge)
            ions_scaled = [<int>(x / acc_frag) for x in marr]
            if iname in theoretical_set:
                theoretical_set_item = <list>PyDict_GetItem(theoretical_set, iname)
                theoretical_set_item.extend(ions_scaled)
            else:
                theoretical_set[iname] = ions_scaled                

            if reshape:
                marr_storage = np.array(marr)
                marr_storage.sort()
                n = marr_storage.size
                marr_storage = marr_storage.reshape((n, 1))

                iname = (ion_type, charge)
                peaks[iname] = marr_storage
            else:
                iname = (ion_type, charge)
                peaks[iname] = sorted(marr)
    return peaks, theoretical_set


def theor_spectrum(peptide, acc_frag, nterm_mass, cterm_mass, types=('b', 'y'), maxcharge=None, reshape=False, **kwargs):
    if not maxcharge:
        maxcharge = 1 + int(ec.charge(peptide, pH=2))
    return ctheor_spectrum(peptide, acc_frag, nterm_mass, cterm_mass, tuple(types), maxcharge, reshape, kwargs)

