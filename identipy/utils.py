import re
from collections import defaultdict
from pyteomics import mass, electrochem as ec, auxiliary as aux, parser
import sys
import numpy as np
from multiprocessing import Queue, Process, cpu_count
from string import punctuation

def decode(func):
    if sys.version_info.major == 3:
        def f(s, *args, **kwargs):
            if isinstance(s, bytes):
                return func(s.decode('ascii'), *args, **kwargs)
            else:
                return func(s, *args, **kwargs)
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

def neutral_masses(spectrum, settings):
    if 'params' in spectrum:
        exp_mass = spectrum['params']['pepmass'][0]
        charge = spectrum['params'].get('charge')
    else:
        ion = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
        charge = int(ion['charge state'])
        exp_mass = ion['selected ion m/z']
    if isinstance(charge, str):
        states = aux._parse_charge(charge)
        if isinstance(states, int):
            states = [states]
    elif charge is None:
        states = list(range(1, 1+settings.getint('search', 'maximum charge')))
    else: states = [charge]
    states = [s for s in states if s <= settings.getint('search', 'maximum charge')]
    states.sort()
    return zip((exp_mass*ch - ch*mass.nist_mass['H+'][0][0]
            for ch in states), states)

@aux.memoize(10)
def import_(name):
    """Import a function by name: module.function or
    module.submodule.function, etc. Return the function object."""

    mod, f = name.rsplit('.', 1)
    return getattr(__import__(mod, fromlist=[f]), f)

def aa_mass(settings):
    if settings.hasoption('misc', 'aa_mass'):
        return settings.get('misc', 'aa_mass')
    aa_mass = mass.std_aa_mass.copy()
    fmods = settings.get('modifications', 'fixed')
    if fmods:
        for mod in re.split(r'[,;]\s*', fmods):
            m, aa = parser._split_label(mod)
            aa_mass[aa] += settings.getfloat('modifications', m)
    vmods = settings.get('modifications', 'variable')
    if vmods:
        mods = [parser._split_label(l) for l in re.split(r',\s*', vmods)]
        for (mod, aa), char in zip(mods, punctuation):
            aa_mass[char] = aa_mass[aa] + settings.getfloat('modifications', mod)

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

def allow_all(*args):
    return True

def get_RT(spectrum):
    if 'rtinseconds' in spectrum:
        return float(spectrum['rtinseconds'])
    return spectrum['scanList']['scan'][0]['scan start time']


#class Config(defaultdict):
#    """
#    A (nested) dict that provides set and get methods to mimic
#    ConfigParser. 
#    """
#    def __init__(self, *args, **kwargs):
#        defaultdict.__init__(self, dict)
#    def getint(self, section, option):
#        return int(self[section][option])
#    def getfloat(self, section, option):
#        return float(self[section][option])
#    def getboolean(self, section, option):
#        val = self[section][option]
#        if isinstance(val, (bool, int, float)):
#            return bool(val)
#        if isinstance(val, str):
#            if val.lower() in {'yes', 'true', '1', 'on'}:
#                return True
#            if val.lower() in {'no', 'false', '0', 'off'}:
#                return False
#        raise TypeError('Cannot convert to boolean: {}'.format(val))
#    def get(self, section, option):
#        return self[section][option]
#    def set(self, section, option, value):
#        self[section][option] = value
#    def has_option(self, section, option):
#        return option in self[section]
#    def has_section(self, section):
#        return section in self
