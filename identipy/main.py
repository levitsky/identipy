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
import scoring, utils
from types import FunctionType
from ConfigParser import RawConfigParser

def candidates_from_arrays(spectrum, settings):
    spectrum = copy(spectrum)
    idx = np.nonzero(spectrum['m/z array'] >=
            settings.getfloat('search', 'product minimum m/z'))
    spectrum['intensity array'] = spectrum['intensity array'][idx]
    spectrum['m/z array'] = spectrum['m/z array'][idx]
    maxpeaks = settings.getint('scoring', 'maximum peaks')
    minpeaks = settings.getint('scoring', 'minimum peaks')
    maxlen = settings.getint('search', 'peptide maximum length')

    dtype = np.dtype([('score', np.float32),
                ('seq', np.str_, maxlen), ('note', np.str_, 1)])
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

    masses, seqs, notes = get_arrays(settings) if not settings.has_option(
            'performance', 'arrays') else settings.get('performance', 'arrays')
    exp_mass = utils.neutral_masses(spectrum, settings)
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

    candidates = []
    candidates_notes = []
    for m, c in exp_mass:
        if c != 1:
            dm_l = acc_l * m / 1.0e6 if rel else acc_l * c
            dm_r = acc_r * m / 1.0e6 if rel else acc_r * c
            start = masses.searchsorted(m - dm_l)
            end = masses.searchsorted(m + dm_r)
            candidates.extend(seqs[start:end])
            candidates_notes.extend(notes[start:end])

    result = []
    for idx, x in enumerate(candidates):
        s = score(spectrum, x, settings)
        if s > 0:
            result.append((s, x, candidates_notes[idx]))
    result.sort(reverse=True)

    if settings.has_option('misc', 'legend'):
        mods = list(zip(settings.get('misc', 'legend'), punctuation))
        res = []
        for score, cand, note in result:
            for (mod, aa), char in mods:
                cand = cand.replace(char, mod + aa)
            if len(cand) > maxlen:
                maxlen = len(cand)
            res.append((score, cand, note))
        result = res
    return np.array(result, dtype=dtype)


def get_arrays(settings):
    db = settings.get('input', 'database')
    hasher = settings.get('misc', 'hash')
    dbhash = hashlib.new(hasher)
    with open(db) as f:
        print "Scanning database contents..."
        for line in f:
            dbhash.update(line)
    dbhash = dbhash.hexdigest()
    print "Done."
    folder = settings.get('performance', 'folder')
    if not os.path.isdir(folder):
        os.makedirs(folder)
    enzyme = settings.get('search', 'enzyme')
    enzyme = parser.expasy_rules.get(enzyme, enzyme)
    mc = settings.getint('search', 'missed cleavages')
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
            func = lambda prot: [(pep, get_note(prot[0], label=prefix))
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
        if not settings.has_option('performance', 'arrays'):
            settings.set('performance', 'arrays', get_arrays(settings))
        candidates = lambda s: candidates_from_arrays(s, settings)
        if processor == 'minimal':
            return lambda s: {'spectrum': s, 'candidates': candidates(s)}
        elif processor == 'e-value':
            condition = settings.get('scoring', 'condition')
            if condition:
                if not isinstance(condition, FunctionType):
                    condition = utils.import_(condition)

            def f(s):
                c = candidates(s)
                if condition:
                    c = np.array([x for x in c if condition(s, x[1], settings)],
                            dtype=c.dtype)
                return {'spectrum': s, 'candidates': c,
                    'e-values': scoring.evalues(c, settings)}
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
    if settings.get('modifications', 'variable'):
        return double_run(fname, settings, varmod_stage1)
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
    new_settings.remove_option('performance', 'arrays')
    print '[double run] stage 2 starting ...'
    return process_file(fname, new_settings)


def varmod_stage1(fname, settings):
    """Take mods, make a function that yields new settings"""
    print 'Running preliminary search (no modifications) ...'
    aa_mass = utils.get_aa_mass(settings)
    mods = settings.get('modifications', 'variable')
    mods = [parser._split_label(l) for l in re.split(r',\s*', mods)]
    mods.sort(key=lambda x: len(x[0]), reverse=True)
    assert all(len(m) == 2 for m in mods), 'unmodified residue given'
    mod_dict = {}
    for mod, aa in mods:
        mod_dict.setdefault(mod, []).append(aa)

    candidates = set()
    settings = copy(settings)
    settings.set('modifications', 'variable', '')
    for s in process_file(fname, settings):
        candidates.update(c[1] for c in s['candidates'])

    n = settings.getint('modifications', 'maximum variable mods')
    seq_iter = chain.from_iterable(
            (parser.tostring(x, False) for x in
                parser.isoforms(seq, variable_mods=mod_dict, format='split')
                if ((len(x[0]) > 2) + (len(x[-1]) > 2) + sum(
                    len(y) > 1 for y in x[1:-1])) <= n)
            for seq in candidates)

    def prepare_seqs():
        for seq in seq_iter:
            for (mod, aa), char in zip(mods, punctuation):
                seq = seq.replace(mod + aa, char)
            yield seq
    maxlen = settings.getint('search', 'peptide maximum length')
    seqs = np.fromiter(prepare_seqs(), dtype=np.dtype((np.str_, maxlen)))
    masses = np.empty(seqs.shape, dtype=np.float32)
    for i in range(seqs.size):
        masses[i] = mass.fast_mass(seqs[i], aa_mass=aa_mass)
    i = masses.argsort()
    masses = masses[i]
    seqs = seqs[i]
    new_settings = copy(settings)
    new_settings.set('misc', 'aa_mass', aa_mass)
    new_settings.set('misc', 'legend', mods)
    new_settings.set('performance', 'arrays', (masses, seqs))
    maxlen = settings.getint('search', 'peptide maximum length')
    return new_settings


@aux.memoize(10)
def settings(fname=None, default_name=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), os.pardir, 'default.cfg')):
    """Read a configuration file and return a :py:class:`RawConfigParser` object.
    """

    raw_config = RawConfigParser(dict_type=dict, allow_no_value=True)
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
