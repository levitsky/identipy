import sys
from pyteomics import mass, electrochem as ec
import numpy as np

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


