import numpy as np
from pyteomics import parser, mass, fasta, auxiliary as aux, mgf, mzml
from itertools import chain
import os
from tempfile import NamedTemporaryFile
import ast
import hashlib
from copy import copy
from string import punctuation
from . import scoring, utils

def candidates_from_arrays(spectrum, settings):
    spectrum = utils.preprocess_spectrum(spectrum, settings)
    
    maxlen = settings.getint('search', 'peptide maximum length')
    dtype = np.dtype([('score', np.float64),
                ('seq', np.str_, maxlen), ('note', np.str_, 1),
                ('charge', np.int8), ('info', np.object_), ('sumI', np.float64)])
    if spectrum is None:
        return np.array([], dtype=dtype)

    masses, seqs, notes = get_arrays(settings)
    exp_mass = utils.neutral_masses(spectrum, settings)
    charge2mass = dict((c, m) for m, c in exp_mass)
    score = utils.import_(settings.get('scoring', 'score'))
    acc_l = settings.getfloat('search', 'precursor accuracy left')
    acc_r = settings.getfloat('search', 'precursor accuracy right')
    unit = settings.get('search', 'precursor accuracy unit')
    rel = utils.relative(unit)
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
        result.append((s.pop('score'), x, candidates_notes[idx], candidates_charges[idx], s, s.pop('sumI'), s.pop('fragmentMT')))
        result[-1][4]['mzdiff'] = {'Th': charge2mass[candidates_charges[idx]] - masses[indexes[idx]]}
        result[-1][4]['mzdiff']['ppm'] = 1e6 * result[-1][4]['mzdiff']['Th'] / masses[indexes[idx]]
    result.sort(reverse=True)

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
    add_decoy = settings.getboolean('input', 'add decoy')
    print "Done."
    folder = settings.get('performance', 'folder')
    if not os.path.isdir(folder):
        os.makedirs(folder)
    enzyme = utils.get_enzyme(settings.get('search', 'enzyme'))
    mc = settings.getint('search', 'number of missed cleavages')
    minlen = settings.getint('search', 'peptide minimum length')
    maxlen = settings.getint('search', 'peptide maximum length')
    aa_mass = utils.get_aa_mass(settings)
    
    index = os.path.join(folder, 'identipy.idx')
    

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
        def get_note(description, label='DECOY_'):
            return 'd' if description.startswith(label) else 't'

        def peps():
            prots = utils.prot_gen(settings)
            aa_mass = utils.get_aa_mass(settings)
            mods = settings.get('modifications', 'variable')
            if mods:
                maxmods = settings.getint('modifications', 'maximum variable mods')
                legend = settings.get('misc', 'legend')
                punct = set(punctuation)
                nmods = [(p, mod[1], mod[2]) for p, mod in legend.iteritems() if p in punct]

                def func(prot):
                    note = get_note(prot[0], label=prefix)
                    out = []
                    for pep in parser.cleave(prot[1], enzyme, mc):
                        if minlen <= len(pep) <= maxlen and parser.fast_valid(pep):
                            for seqm in utils.custom_isoforms(pep, variable_mods=nmods, maxmods=maxmods):
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
        utils.set_mod_dict(settings)

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



