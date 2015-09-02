import numpy as np
from pyteomics import parser, mass, fasta, auxiliary as aux, mgf, mzml
from itertools import chain
import re
import os
from tempfile import NamedTemporaryFile
import ast
import hashlib
from copy import copy
from string import punctuation
from . import scoring, utils
from utils import CustomRawConfigParser, get_enzyme
# from ConfigParser import RawConfigParser
import operator as op

try:
    from pyteomics import cmass
except ImportError:
    cmass = mass

def candidates_from_arrays(spectrum, settings):
    spectrum = copy(spectrum)
    idx = np.nonzero(spectrum['m/z array'] >=
            settings.getfloat('search', 'product minimum m/z'))
    spectrum['intensity array'] = spectrum['intensity array'][idx]
    spectrum['m/z array'] = spectrum['m/z array'][idx]
    maxpeaks = settings.getint('scoring', 'maximum peaks')
    minpeaks = settings.getint('scoring', 'minimum peaks')
    maxlen = settings.getint('search', 'peptide maximum length')

    dtype = np.dtype([('score', np.float64),
                ('seq', np.str_, maxlen), ('note', np.str_, 1),
                ('charge', np.int8), ('info', np.object_), ('sumI', np.float64)])
    if maxpeaks and minpeaks > maxpeaks:
        raise ValueError('minpeaks > maxpeaks: {} and {}'.format(
            minpeaks, maxpeaks))
    if maxpeaks and spectrum['intensity array'].size > maxpeaks:
        i = np.argsort(spectrum['intensity array'])[-maxpeaks:]
        j = np.argsort(spectrum['m/z array'][i])
        spectrum['intensity array'] = spectrum['intensity array'][i][j]
        spectrum['m/z array'] = spectrum['m/z array'][i][j]
    if minpeaks and spectrum['intensity array'].size < minpeaks:
        return np.array([], dtype=dtype)
    dynrange = settings.getfloat('scoring', 'dynamic range')
    if dynrange:
        i = spectrum['intensity array'] > spectrum['intensity array'].max(
                ) / dynrange
        spectrum['intensity array'] = spectrum['intensity array'][i]
        spectrum['m/z array'] = spectrum['m/z array'][i]
    if minpeaks and spectrum['intensity array'].size < minpeaks:
        return np.array([], dtype=dtype)

    masses, seqs, notes = get_arrays(settings)
    exp_mass = utils.neutral_masses(spectrum, settings)
    charge2mass = dict((c, m) for m, c in exp_mass)
    score = utils.import_(settings.get('scoring', 'score'))
    acc_l = settings.getfloat('search', 'precursor accuracy left')
    acc_r = settings.getfloat('search', 'precursor accuracy right')
    unit = settings.get('search', 'precursor accuracy unit')
    if unit == 'ppm':
        rel = True
    elif unit in {'Th', 'Da', 'amu'}:
        rel = False
    else:
        raise ValueError('Unrecognized precursor accuracy unit: ' + unit)
    aa_mass = utils.get_aa_mass(settings)

    candidates = []
    candidates_notes = []
    candidates_charges = []
    candidates_info = []
    indexes = []
    for m, c in exp_mass:
        dm_l = acc_l * m / 1.0e6 if rel else acc_l * c
        dm_r = acc_r * m / 1.0e6 if rel else acc_r * c
        start = masses.searchsorted(m - dm_l)
        end = masses.searchsorted(m + dm_r)
        candidates.extend(seqs[start:end])
        candidates_notes.extend(notes[start:end])
        candidates_charges.extend([c] * (end-start))
        indexes.extend(range(start, end))

    result = []
    for idx, x in enumerate(candidates):
        s = score(spectrum, x, candidates_charges[idx], settings)
        result.append((s.pop('score'), x, candidates_notes[idx], candidates_charges[idx], s, s.pop('sumI')))
        result[-1][4]['mzdiff'] = {'Th': charge2mass[candidates_charges[idx]] - masses[indexes[idx]]}
        result[-1][4]['mzdiff']['ppm'] = 1e6 * result[-1][4]['mzdiff']['Th'] / masses[indexes[idx]]
    result.sort(reverse=True)

    if settings.has_option('misc', 'legend'):
        leg = settings.get('misc', 'legend')
        res = []
        for score, cand, note, charge, info, sumI in result:
            for char in punctuation:
                if char in leg:
                    cand = cand.replace(char, ''.join(leg[char]))
            if len(cand) > maxlen:
                maxlen = len(cand)
            res.append((score, cand, note, charge, info, sumI))
        result = res
    return np.array(result, dtype=dtype)


def get_arrays(settings):
    if settings.has_option('performance', 'arrays'):
        return settings.get('performance', 'arrays')
    db = settings.get('input', 'database')
    hasher = settings.get('misc', 'hash')
    dbhash = hashlib.new(hasher)
    with open(db) as f:
        print "Scanning database contents ..."
        for line in f:
            dbhash.update(line)
    dbhash = dbhash.hexdigest()
    print "Done."
    folder = settings.get('performance', 'folder')
    if not os.path.isdir(folder):
        os.makedirs(folder)
    enzyme = settings.get('search', 'enzyme')
    enzyme = get_enzyme(enzyme)
    mc = settings.getint('search', 'number of missed cleavages')
    minlen = settings.getint('search', 'peptide minimum length')
    maxlen = settings.getint('search', 'peptide maximum length')
    aa_mass = utils.get_aa_mass(settings)
    add_decoy = settings.getboolean('input', 'add decoy')
    index = os.path.join(folder, 'identipy.idx')
    prefix = settings.get('input', 'decoy prefix')

    profile = (dbhash, add_decoy, enzyme, mc, minlen, maxlen, aa_mass)
    arr_name = None
    if os.path.isfile(index):
        with open(index) as ifile:
            for line in ifile:
                t, name = line.strip().split('\t')
                if ast.literal_eval(t) == profile:
                    arr_name = name
                    break

    if arr_name is None:

        print 'Generating peptide arrays ...'
        def get_note(protein_description, label='DECOY_'):
            return 'd' if protein_description.startswith(label) else 't'

        def peps():
            if not add_decoy:
                prots = fasta.read(db)
            else:
                mode = settings.get('input', 'decoy method')
                prots = fasta.decoy_db(db, mode=mode, prefix=prefix)

            aa_mass = utils.get_aa_mass(settings)
            mods = settings.get('modifications', 'variable')
            if mods:
                maxmods = settings.getint('modifications', 'maximum variable mods')
                legend = settings.get('misc', 'legend')
                def func(prot):
                    note = get_note(prot[0], label=prefix)
                    out = []
                    for pep in parser.cleave(prot[1], enzyme, mc):
                        if minlen <= len(pep) <= maxlen and parser.fast_valid(pep):
                            for seq in parser.isoforms(pep, variable_mods=mods, maxmods=maxmods):
                                seqm = seq
                                for res, char in legend.iteritems():
                                    seqm = seqm.replace(res, char)
                                out.append((seqm, note))
                    return out
            else:
                def func(prot):
                    note = get_note(prot[0], label=prefix)
                    return [(pep, note)
                        for pep in parser.cleave(prot[1], enzyme, mc)
                        if minlen <= len(pep) <= maxlen and parser.fast_valid(pep)]

            n = settings.getint('performance', 'processes')
            return chain.from_iterable(
                    utils.multimap(n, func, prots))

        seqs_and_notes = np.fromiter(peps(), dtype=np.dtype(
            [('seq', np.str_, maxlen), ('note', np.str_, 1)]))

        seqs, idx = np.unique(seqs_and_notes['seq'], return_index=True)
        notes = seqs_and_notes['note'][idx]
        del seqs_and_notes
        masses = np.empty(seqs.shape, dtype=np.float32)
        for i in np.arange(seqs.size):
            masses[i] = mass.fast_mass(seqs[i], aa_mass=aa_mass)
        idx = np.argsort(masses)
        masses = masses[idx]
        seqs = seqs[idx]
        notes = notes[idx]
        with NamedTemporaryFile(suffix='.npz', prefix='identipy_', dir=folder,
                delete=False) as outfile:
            np.savez_compressed(outfile, masses=masses, seqs=seqs, notes=notes)
            name = outfile.name
        with open(index, 'a') as ifile:
            ifile.write(str(profile) + '\t' + name + '\n')
        print "Arrays saved."

    else:
        npz = np.load(arr_name)
        seqs = npz['seqs']
        masses = npz['masses']
        notes = npz['notes']
    return masses, seqs, notes


def spectrum_processor(settings):
    processor = settings.get('misc', 'spectrum processor')
    if '.' in processor:
        return utils.import_(processor)(settings)

    mode = settings.get('performance', 'pre-calculation')
    if mode == 'some':  # work with numpy arrays
        settings = copy(settings)
        
        mods = settings.get('modifications', 'variable').strip()
        mod_dict = {}
        if mods:
            legend = {}
            mods = [parser._split_label(l) for l in re.split(r',\s*', mods)]
            mods.sort(key=lambda x: len(x[0]), reverse=True)
            for mod, char in zip(mods, punctuation):
                legend[''.join(mod)] = char
                legend[char] = mod
            assert all(len(m) == 2 for m in mods), 'unmodified residue given'
            for mod, aa in mods:
                mod_dict.setdefault(mod, []).append(aa)
            settings.set('misc', 'legend', legend)
        settings.set('modifications', 'variable', mod_dict)

        if not settings.has_option('performance', 'arrays'):
            settings.set('performance', 'arrays', get_arrays(settings))
        candidates = lambda s: candidates_from_arrays(s, settings)
        if processor == 'minimal':
            return lambda s: {'spectrum': s, 'candidates': candidates(s)}

        elif processor.startswith('e-value'):
            condition = settings.get('scoring', 'condition')
            if condition:
                if not callable(condition):
                    condition = utils.import_(condition)

            def f(s):
                c = candidates(s)
                if condition:
                    c = np.array([tuple(x) for x in c if condition(s, x[1], settings)],
                            dtype=c.dtype)
                if processor == 'e-value':
                    return {'spectrum': s, 'candidates': c,
                        'e-values': scoring.evalues(c, settings)}
                elif processor == 'e-value2':
                    if not hasattr(spectrum_processor, 'cache'):
                        spectrum_processor.cache = {}
                    for cand in c[:]:
                        if cand[0] > 0 and cand[2] == 'd':
                            for m, ch in utils.neutral_masses(s, settings):
                                spectrum_processor.cache.setdefault(
                                        (len(cand[1]), ch), []).append(cand[0])
                    return {'spectrum': s, 'candidates': c}
#                       'e-values': scoring.evalues2(s, c, settings)}
            return f
    else:
        raise NotImplementedError('Unsupported pre-calculation mode')


def process_spectra(f, settings):
    # prepare the function
    func = spectrum_processor(settings)
    print 'Running the search ...'
    # decide on multiprocessing
    n = settings.getint('performance', 'processes')
    return utils.multimap(n, func, f)

def process_file(fname, settings):
    stage1 = settings.get('misc', 'first stage')
    if stage1:
        return double_run(fname, settings, utils.import_(stage1))
    else:
        ftype = fname.rsplit('.', 1)[-1].lower()
        if ftype == 'mgf':
            spectra = mgf.read(fname)
        elif ftype == 'mzml':
            spectra = (x for x in mzml.read(fname) if x['ms level'] > 1)
        else:
            raise ValueError('Unrecognized file type: {}'.format(ftype))
        return process_spectra(spectra, settings)

def double_run(fname, settings, stage1):
    print '[double run] stage 1 starting ...'
    new_settings = stage1(fname, settings)
    print '[double run] stage 2 starting ...'
    return process_file(fname, new_settings)


def settings(fname=None, default_name=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'default.cfg')):
    """Read a configuration file and return a :py:class:`RawConfigParser` object.
    """

    raw_config = CustomRawConfigParser(dict_type=dict, allow_no_value=True)
    if default_name:
        raw_config.read(default_name)
    if fname:
        raw_config.read(fname)
    return raw_config
#   config = utils.Config()
#   conv = {True: utils.import_, False: ast.literal_eval}
#   for section in raw_config.sections():
#       for option in raw_config.options(section):
#           val = raw_config.get(section, option)
#           config[section][option] = conv[
#                   raw_config.has_option('import', section + ', ' + option)](
#                           raw_config.get(section, option)) if val else None
#   return config
