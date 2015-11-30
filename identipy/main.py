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
from collections import defaultdict
import operator as op
from bisect import bisect
from . import scoring, utils

try:
    from pyteomics import cmass
except ImportError:
    cmass = mass

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

def prepare_peptide_processor(fname, settings):
    global spectra
    global nmasses
    global charges
    global idx
    global t2i
    spectra = []
    nmasses = []
    charges = []
    idx = []
    t2i = {}
    print 'Reading spectra ...'
    for spec in iterate_spectra(fname):
        ps = utils.preprocess_spectrum(spec, settings)
        if ps is not None:
            spectra.append(ps)
    for i, s in enumerate(spectra):
        t2i[utils.get_title(s)] = i
        for m, c in utils.neutral_masses(s, settings):
            charges.append(c)
            nmasses.append(m)
            s.setdefault('nm', []).append(m)
            s.setdefault('ch', []).append(c)
            idx.append(i)
    print len(spectra), 'spectra pass quality criteria.'
    i = np.argsort(nmasses)
    nmasses = np.array(nmasses)[i]
    charges = np.array(charges)[i]
    idx = np.array(idx)[i]
    spectra = np.array(spectra)

    utils.set_mod_dict(settings)

    mods = settings.get('modifications', 'variable')
    maxmods = settings.getint('modifications', 'maximum variable mods')
    leg = settings.get('misc', 'legend')
    punct = set(punctuation)
    nmods = [(p, mod[1], mod[2]) for p, mod in leg.iteritems() if p in punct]

    aa_mass = utils.get_aa_mass(settings)
    score = utils.import_(settings.get('scoring', 'score'))
    acc_l = settings.getfloat('search', 'precursor accuracy left')
    acc_r = settings.getfloat('search', 'precursor accuracy right')
    acc_frag = settings.getfloat('search', 'product accuracy')
    unit = settings.get('search', 'precursor accuracy unit')
    rel = utils.relative(unit)
    return {'rel': rel, 'aa_mass': aa_mass, 'acc_l': acc_l, 'acc_r': acc_r, 'acc_frag': acc_frag,
            'unit': unit, 'nmods': nmods, 'maxmods': maxmods, 'settings': settings}

def peptide_processor_iter_isoforms(peptide, **kwargs):
    nmods, maxmods = op.itemgetter('nmods', 'maxmods')(kwargs)
    out = []
    if nmods and maxmods:
        for form in utils.custom_isoforms(peptide, variable_mods=nmods, maxmods=maxmods):
            out.append(peptide_processor(form, **kwargs))
    else:
        out.append(peptide_processor(peptide, **kwargs))
    return out


def peptide_processor(peptide, **kwargs):
    seqm = peptide
    m = cmass.fast_mass(seqm, aa_mass=kwargs['aa_mass'])
    rel = kwargs['rel']
    acc_l = kwargs['acc_l']
    acc_r = kwargs['acc_r']
    settings = kwargs['settings']
    dm_l = acc_l * m / 1.0e6 if rel else acc_l# * c FIXME
    dm_r = acc_r * m / 1.0e6 if rel else acc_r# * c FIXME
    start = nmasses.searchsorted(m - dm_l)
    end = nmasses.searchsorted(m + dm_r)
    if start == end: return None
    cand_idx = idx[start:end]
    cand_spectra = spectra[cand_idx]
    if settings.has_option('scoring', 'condition'):
        cond = settings.get('scoring', 'condition')
    else:
        cond = None
    if isinstance(cond, str) and cond.strip():
        cond = utils.import_(cond)
    if cond:
        cand_spectra = [c for c in cand_spectra if cond(c, seqm, settings)]
    theor = utils.theor_spectrum(seqm, maxcharge=1, aa_mass=kwargs['aa_mass'])
    results = [scoring._hyperscore(copy(s), theor, kwargs['acc_frag']) for s in cand_spectra] # FIXME (use score from settings?)
    
    results = [(x.pop('score'), utils.get_title(s), x, m) for x, s in zip(results, cand_spectra)]

    results.sort(reverse=True)
    # results = np.array(results, dtype=[('score', np.float32), ('title', np.str_, 30), ('spectrum', np.object_), ('info', np.object_)])
    return peptide, results


def iterate_spectra(fname):
    ftype = fname.rsplit('.', 1)[-1].lower()
    if ftype == 'mgf':
        with mgf.read(fname) as f:
            for x in f:
                yield x
    elif ftype == 'mzml':
        with mzml.read(fname) as f:
            for x in f:
                if x['ms level'] > 1:
                    yield x
    else:
        raise ValueError('Unrecognized file type: {}'.format(ftype))

def process_peptides(fname, settings):
    # global spec_scores
    # global spec_top_scores
    # global spec_top_seqs

    spec_results = defaultdict(dict)
    spec_scores = {}
    spec_top_scores = {}
    spec_top_seqs = {}
    peps = utils.peptide_gen(settings)
    kwargs = prepare_peptide_processor(fname, settings)
    func = peptide_processor_iter_isoforms
    nc = settings.getint('output', 'candidates') or None
    print 'Running the search ...'
    n = settings.getint('performance', 'processes')
    leg = {}
    if settings.has_option('misc', 'legend'):
        leg = settings.get('misc', 'legend')
    for y in utils.multimap(n, func, peps, **kwargs):
        for x in y:
            if x is not None:
                peptide, result = x
                for score, spec_t, info, m in result:
                    spec = spectra[t2i[spec_t]]
                    score = float(score)
                    info['pep_nm'] = m
                    spec_results[spec_t]['spectrum'] = spec
#               spec_results[spec_t].setdefault('scores', []).append(score) FIXME write histogram
                    spec_results[spec_t].setdefault('sequences', [])
                    spec_results[spec_t].setdefault('top_scores', [])
                    spec_results[spec_t].setdefault('info', [])
                    # spec_scores.setdefault(spec_t, []).append(score)
                    # spec_top_scores.setdefault(spec_t, [])
                    # spec_top_seqs.setdefault(spec_t, [])
                    top_scores = spec_results[spec_t]['top_scores']
                    top_seqs = spec_results[spec_t]['sequences']
                    top_info = spec_results[spec_t]['info']
                    i = bisect(top_scores, -score)
                    if nc is None or i < nc:
                        top_scores.insert(i, -score)
                        top_seqs.insert(i, peptide)
                        top_info.insert(i, info)
                        if nc is not None and len(top_scores) > nc:
                            top_scores.pop()
                            top_seqs.pop()
                            top_info.pop()
    maxlen = settings.getint('search', 'peptide maximum length')
    dtype = np.dtype([('score', np.float64),
        ('seq', np.str_, maxlen), ('note', np.str_, 1),
        ('charge', np.int8), ('info', np.object_), ('sumI', np.float64), ('fragmentMT', np.float64)])

    for spec_name, val in spec_results.iteritems():
        s = val['spectrum']
        c = []
        evalues = []
        for idx, score in enumerate(val['top_scores']):
            mseq = val['sequences'][idx]
            seq = mseq
            for x in set(mseq).intersection(punctuation):
                seq = seq.replace(x, leg[x][1])
            pnm = val['info'][idx]['pep_nm']
            nidx = min(range(len(s['nm'])), key=lambda i: abs(s['nm'][i]-pnm))
            c.append((-score, mseq, 't' if seq in utils.seen_target else 'd', s['ch'][nidx], val['info'][idx], val['info'][idx].pop('sumI'), val['info'][idx].pop('fragmentMT')))
            c[-1][4]['mzdiff'] = {'Th': s['nm'][nidx] - pnm}
            c[-1][4]['mzdiff']['ppm'] = 1e6 * c[-1][4]['mzdiff']['Th'] / pnm
            evalues.append(-1./score if -score else 1e6)
        c = np.array(c, dtype=dtype)
        yield {'spectrum': s, 'candidates': c, 'e-values': evalues}

    # return spec_results
    # return (func(p) for p in peps)
    # return utils.multimap(n, func, peps)

def process_file(fname, settings):
    
    stage1 = settings.get('misc', 'first stage')
    if stage1:
        return double_run(fname, settings, utils.import_(stage1))
    else:
        iterate = settings.get('misc', 'iterate')
        ftype = fname.rsplit('.', 1)[-1].lower()
        if iterate == 'spectra':
            spectra = iterate_spectra(fname)
            return process_spectra(spectra, settings)
        elif iterate == 'peptides':
            utils.seen_target.clear()
            utils.seen_decoy.clear()
            return process_peptides(fname, settings)
        else:
            raise ValueError('iterate must be "spectra" or "peptides"')

def double_run(fname, settings, stage1):
    print '[double run] stage 1 starting ...'
    new_settings = stage1(fname, settings)
    print '[double run] stage 2 starting ...'
    return process_file(fname, new_settings)


def settings(fname=None, default_name=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'default.cfg')):
    """Read a configuration file and return a :py:class:`RawConfigParser` object.
    """

    raw_config = utils.CustomRawConfigParser(dict_type=dict, allow_no_value=True)
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
