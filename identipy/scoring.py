from .utils import neutral_masses
from .analysis import theor_spectrum
import numpy as np
acc = 0.02

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


