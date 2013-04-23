from .utils import neutral_masses, theor_spectrum, aa_mass
from scipy.spatial import cKDTree
import numpy as np
from scipy.stats import scoreatpercentile
from math import factorial

def simple_score(spectrum, peptide, settings):
    acc = settings.getfloat('search', 'product accuracy')
    charge = max(c for m, c in neutral_masses(spectrum))
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
    int_array = int_array / int_array.max()
    acc = settings.getfloat('search', 'product accuracy')
    charge = max(c for m, c in neutral_masses(spectrum))
    theor = theor_spectrum(peptide, maxcharge=charge, aa_mass=aa_mass(settings))
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

def bin_size_sorted(X):
# from http://mail.scipy.org/pipermail/scipy-user/2009-March/020194.html
    """Calculates the Freedman-Diaconis bin size for a sorted data set.

    Paramters:
    ----------
        X:  1D data set
    Returns:
    --------
        h:  F-D bin size
    """
#   X = np.sort(X)

    upperQuartile = scoreatpercentile(X, 75)
    lowerQuartile = scoreatpercentile(X, 25)
    IQR = upperQuartile - lowerQuartile

    h = 2. * IQR / len(X)**(1./3.)
    return h
