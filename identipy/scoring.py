from .utils import neutral_masses, theor_spectrum, get_aa_mass
from scipy.spatial import cKDTree
import numpy as np
from math import factorial
from copy import copy


def simple_score(spectrum, peptide, settings):
    acc = settings.getfloat('search', 'product accuracy')
    int_array = spectrum['intensity array']
    int_array = int_array / int_array.max() * 100
    charge = max(c for _, c in neutral_masses(spectrum, settings))
    theor = theor_spectrum(peptide, maxcharge=charge, aa_mass=get_aa_mass(settings))
    fragments = np.concatenate(theor.values())
    if '__KDTree' not in spectrum:
        spectrum['__KDTree'] = cKDTree(spectrum['m/z array'].reshape(
            (spectrum['m/z array'].size, 1)))
    dist, ind = spectrum['__KDTree'].query(fragments.reshape((fragments.size, 1)),
            distance_upper_bound=acc)
    mask = dist != np.inf
    if mask.size < settings.getint('scoring', 'minimum matched'):
        return -1
    return spectrum['intensity array'][ind[mask]].sum()


def get_fragment_mass_tol(spectrum, peptide, settings):
    """A function for obtaining optimal fragment mass tolerance, dynamic range"""
    acc = 1.0  # maximum product accuracy
    spectrum = copy(spectrum)
    idx = np.nonzero(spectrum['m/z array'] >= 150)
    spectrum['intensity array'] = spectrum['intensity array'][idx]
    spectrum['m/z array'] = spectrum['m/z array'][idx]
    int_array = spectrum['intensity array']
    int_array = int_array / int_array.max() * 100
    charge = max(1, max(c for _, c in neutral_masses(spectrum, settings)) - 1)
    theor = theor_spectrum(peptide, maxcharge=charge, aa_mass=get_aa_mass(settings))
    if '__KDTree' not in spectrum:
        spectrum['__KDTree'] = cKDTree(spectrum['m/z array'].reshape(
            (spectrum['m/z array'].size, 1)))

    dist_total, int_array_total = np.array([]), np.array([])
    for fragments in theor.values():
        n = fragments.size
        dist, ind = spectrum['__KDTree'].query(fragments.reshape((n, 1)),
            distance_upper_bound=acc)
        mask = (dist != np.inf)
        int_array_total = np.append(int_array_total, int_array[ind[mask]])
        dist_total = np.append(dist_total, dist[dist != np.inf])

    new_params = {}
    if dist_total.size:
        new_params['fmt'] = dist_total
    else:
        new_params['fmt'] = None
    if int_array_total.size:
        new_params['dynamic range'] = int_array_total
    else:
        new_params['dynamic range'] = []
    if new_params['dynamic range'].size:
        new_params['maximum peaks'] = int_array_total
    return new_params


def morpheusscore(spectrum, peptide, settings):
    """A simple implementation of Morpheus's score."""
    int_array = spectrum['intensity array']
    int_array = int_array / int_array.max() * 100
    acc = settings.getfloat('search', 'product accuracy')
    charge = max(1, max(c for _, c in neutral_masses(spectrum, settings)) - 1)
    theor = theor_spectrum(peptide, maxcharge=charge, aa_mass=get_aa_mass(settings))
    score = 0
    total_matched = 0
    if '__KDTree' not in spectrum:
        spectrum['__KDTree'] = cKDTree(spectrum['m/z array'].reshape(
            (spectrum['m/z array'].size, 1)))

    for fragments in theor.values():
        n = fragments.size
        dist, ind = spectrum['__KDTree'].query(fragments.reshape((n, 1)),
            distance_upper_bound=acc)
        mask = (dist != np.inf)
        total_matched += fragments.size
        score += int_array[ind[mask]].sum()
    if total_matched < settings.getint('scoring', 'minimum matched'):
        return -1
    return total_matched + score / int_array.sum()


def hyperscore(spectrum, peptide, charge, settings):
    """A simple implementation of X!Tandem's Hyperscore."""
    int_array = spectrum['intensity array']
    mz_array = spectrum['m/z array']
    int_array = int_array / int_array.max() * 100
    acc = settings.getfloat('search', 'product accuracy')
#   charge = max(1, max(c for m, c in neutral_masses(spectrum, settings)) - 1)
    theor = theor_spectrum(peptide, maxcharge=1,#max(1, charge-1),
            aa_mass=get_aa_mass(settings))
    score = 0
    mult = []
    total_matched = 0
    if '__KDTree' not in spectrum:
        spectrum['__KDTree'] = cKDTree(mz_array.reshape((mz_array.size, 1)))

    for fragments in theor.values():
        n = fragments.size
        dist, ind = spectrum['__KDTree'].query(fragments.reshape((n, 1)),
            distance_upper_bound=acc)
        mask = (dist != np.inf)
        nmatched = mask.sum()
        mult.append(factorial(nmatched))
        total_matched += nmatched
        score += int_array[ind[mask]].sum()
    if total_matched < settings.getint('scoring', 'minimum matched'):
        return -1
    for m in mult:
        score *= m

    return score


def survival_hist(scores):
    X_axis = Y_axis = None
    calib_coeff = (-0.18, 3.5)
    if scores.shape[0] > 20:
        hyperscore_h, _ = np.histogram(scores, bins=np.arange(0, round(scores[0]) + 1.5))
        j = hyperscore_h.size - 1
        if j > 20:
            survival_h = hyperscore_h.sum() - np.hstack(([0], hyperscore_h[:-1].cumsum()))
            surv_left = survival_h[0] / 5.
            decr = 0
            while j > 0:
                if survival_h[j] == survival_h[j - 1] and survival_h[j] <= surv_left:
                    decr = survival_h[j]
                    j -= 1
                    while (survival_h[j] == decr and survival_h[j] <= surv_left):
                        survival_h[j] -= decr
                        j -= 1
                else:
                    survival_h[j] -= decr
                    j -= 1
            survival_h[0] -= decr
            max_surv = survival_h[0] * 0.5 + 1.
            min_surv = 10
            proper_surv = (min_surv <= survival_h) * (survival_h <= max_surv)
            if proper_surv.sum() > 2:
                X_axis = proper_surv.nonzero()[0]
                Y_axis = np.log10(survival_h[proper_surv])
                calib_coeff = np.polyfit(X_axis, Y_axis, 1)

    return (X_axis, Y_axis), calib_coeff


def evalues(candidates, settings):
    n = settings.getint('scoring', 'e-values for candidates')
    scores = 4. * np.log10(candidates['score'])
    calib_coeff = survival_hist(scores)[1]
    return 10 ** (scores[:n] * calib_coeff[0] + calib_coeff[1])
