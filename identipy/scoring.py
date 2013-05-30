from .utils import neutral_masses, theor_spectrum, aa_mass, _get_types
from scipy.spatial import cKDTree
import numpy as np
from scipy.stats import scoreatpercentile
from math import factorial

def simple_score(spectrum, peptide, settings):
    acc = settings.getfloat('search', 'product accuracy')
    charge = max(c for m, c in neutral_masses(spectrum, settings))
    theor = theor_spectrum(peptide, maxcharge=charge, aa_mass=aa_mass(settings))
    fragments = np.concatenate(theor.values())
    dist, ind = spectrum['__KDTree'].query(fragments.reshape((fragments.size, 1)),
            distance_upper_bound = acc)
    mask = dist != np.inf
    if mask.size < settings.getint('scoring', 'minimum matched'):
        return -1
    return spectrum['intensity array'][ind[mask]].sum()

def hyperscore(spectrum, peptide, settings):
    """A simple implementation of X!Tandem's Hyperscore."""
    int_array = spectrum['intensity array']
    int_array = 100. * int_array / int_array.max()
    acc = settings.getfloat('search', 'product accuracy')
    charge = max(c for m, c in neutral_masses(spectrum, settings))
    typestr = settings.get('scoring', 'ion types')
    types = _get_types(typestr)
    theor = theor_spectrum(peptide, types=types, 
            maxcharge=charge, aa_mass=aa_mass(settings))
    score = 0
    mult = []
    total_matched = 0
    if '__KDTree' not in spectrum:
        spectrum['__KDTree'] = cKDTree(spectrum['m/z array'].reshape(
            (spectrum['m/z array'].size, 1)))

    for fragments in theor.values():
        n = fragments.size
        dist, ind = spectrum['__KDTree'].query(fragments.reshape((n, 1)),
            distance_upper_bound=acc, eps=acc)
        mask = (dist != np.inf)
        mult.append(factorial(mask.sum()))
        total_matched += fragments.size
        score += int_array[ind[mask]].sum()
    if total_matched < settings.getint('scoring', 'minimum matched'):
        return -1
    if score:
        for m in mult: score *= m
    return score

#def bin_size_sorted(X):
## from http://mail.scipy.org/pipermail/scipy-user/2009-March/020194.html
#    """Calculates the Freedman-Diaconis bin size for a sorted data set.
#
#    Paramters:
#    ----------
#        X:  1D data set
#    Returns:
#    --------
#        h:  F-D bin size
#    """
##   X = np.sort(X)
#
#    upperQuartile = scoreatpercentile(X, 75)
#    lowerQuartile = scoreatpercentile(X, 25)
#    IQR = upperQuartile - lowerQuartile
#
#    h = 2. * IQR / len(X)**(1./3.)
#    return h

def evalues(candidates, settings):
    n = settings.getint('scoring', 'e-values for candidates')
    scores = 4. * np.log10(np.array([x[0] for x in candidates]))
    if len(candidates) < 20:
        calib_coeff = (-0.18, 3.5, 1, 0)
    else:
#       scores.sort() # they are already sorted descending; doesn't matter for histograms
        hyperscore_h, hyperscore_bins = np.histogram(scores, bins=np.arange(0, round(scores[0]) + 1.5))
#       survival_h = [sum(hyperscore_h[x:]) for x in range(len(hyperscore_h))]
        survival_h = hyperscore_h.sum() - np.hstack(([0], hyperscore_h[:-1].cumsum()))
        surv_left = survival_h[0] / 5.
        decr = 0
        j = len(survival_h) - 1
        while j > 0:
            if survival_h[j] == survival_h[j-1] and survival_h[j] <= surv_left:
                decr = survival_h[j]
                j -= 1
                while (survival_h[j] == decr and survival_h[j] <= surv_left):
                    survival_h[j] -= decr
                    j -= 1
            else:
                survival_h[j] -= decr
                j -= 1
        survival_h[0] -= decr

        if len(survival_h) > 20:
            max_surv = survival_h[0] / 2. + 1.
            min_surv = 10
            proper_surv = (min_surv <= survival_h ) * (survival_h <= max_surv)
            if proper_surv.sum() < 2:
                calib_coeff = (-0.18, 3.5, 1, 0)
            else:
                X_axis = proper_surv.nonzero()[0]
                Y_axis = np.log10(survival_h[proper_surv])
#               X_axis, Y_axis = zip(*[(idx, log10(x)) for idx, x in enumerate(survival_h) if min_surv <= x <= max_surv])
#               calib_coeff = auxiliary.linear_regression(X_axis, Y_axis)
                calib_coeff = np.polyfit(X_axis, Y_axis, 1)
        else:
            calib_coeff = (-0.18, 3.5, 1, 0)
#   return [10**(calib_coeff[0] * x + calib_coeff[1]) for x in scores[-n:]]
    return 10 ** (scores[:n] * calib_coeff[0] + calib_coeff[1])
