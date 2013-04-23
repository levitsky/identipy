import re
from pyteomics import mass, electrochem as ec, auxiliary as aux, parser
import sys
import numpy as np
from multiprocessing import Queue, Process, cpu_count
try:
    from configparser import RawConfigParser
except ImportError:
    from ConfigParser import RawConfigParser

def decode(func):
    if sys.version_info.major == 3:
        def f(s, *args, **kwargs):
            if isinstance(s, bytes):
                return func(s.decode('ascii'), *args, **kwargs)
        return f
    return func

@decode
def theor_spectrum(peptide, types=('y', 'b'), maxcharge=None, **kwargs):
    peaks = {}
    if not maxcharge:
        maxcharge = 1 + int(ec.charge(peptide, pH=2))
    for ion_type in types:
        ms = []
        for i in range(1, len(peptide)-1):
            for charge in range(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    ms.append(mass.fast_mass(
                        str(peptide)[:i], ion_type=ion_type, charge=charge,
                        **kwargs))
                else:
                    ms.append(mass.fast_mass(
                        str(peptide)[i:], ion_type=ion_type, charge=charge,
                        **kwargs))
        marr = np.array(ms)
        marr.sort()
        peaks[ion_type] = marr
    return peaks

def neutral_masses(spectrum):
    if 'params' in spectrum:
        exp_mass = spectrum['params']['pepmass'][0]
        charge = spectrum['params'].get('charge')
    else:
        exp_mass = spectrum['base peak']
        charge = spectrum['charge']
    if isinstance(charge, str):
        states = aux._parse_charge(charge)
        if isinstance(states, int):
            states = [states]
    else: states = [charge]
    states.sort()
    return zip((exp_mass*ch - ch*mass.nist_mass['H+'][0][0]
            for ch in states), states)
@aux.memoize(10)
def import_function(name):
    """Import a function by name: module.function or
    module.submodule.function, etc. Return the function object."""

    mod, f = name.rsplit('.', 1)
    return getattr(__import__(mod, fromlist=[f]), f)

@aux.memoize(10)
def aa_mass(settings):
    aa_mass = mass.std_aa_mass.copy()
    fmods = settings.get('modifications', 'fixed')
    if fmods:
        for mod in re.split(r'[,;]\s*', fmods):
            m, aa = parser._split_label(mod)
            aa_mass[aa] += settings.getfloat('modifications', m)
    return aa_mass

def multimap(n, func, it):
    if n == 0:
        try:
            n = cpu_count()
        except NotImplementedError:
            n = 1
    if n == 1:
        for s in it:
            yield func(s)
    else:
        def worker(qin, qout):
            for spectrum in iter(qin.get, None):
                result = func(spectrum)
                qout.put(result)
        qin = Queue()
        qout = Queue()
        count = 0
        for _ in range(n):
            Process(target=worker, args=(qin, qout)).start()
        for s in it:
            qin.put(s)
            count += 1
        for _ in range(n):
            qin.put(None)
        while count:
            yield qout.get()
            count -= 1

class Config(RawConfigParser):
    """
    A ConfigParser with get methods short-circuited to the dict access.
    `optionxform` is not called. `KeyError` will be raised instead of
    `NoSectionError` and `NoOptionError`.
    """
    def __init__(self, *args, **kwargs):
       RawConfigParser.__init__(self, *args, **kwargs)
       self._update_hash()
    def getint(self, section, option):
        return int(self._sections[section][option])
    def getfloat(self, section, option):
        return float(self._sections[section][option])
    def get(self, section, option, **kwargs):
        return self._sections[section][option]
    def set(self, section, option, value):
        RawConfigParser.set(self, section, option, value)
        self._update_hash()
    def __hash__(self):
        return self._hash
    def _get_hash(self):
        return hash(
                frozenset((s, frozenset(opts.items()))
                    for s, opts in self._sections.items()))
    def _update_hash(self):
        self._hash = self._get_hash()
