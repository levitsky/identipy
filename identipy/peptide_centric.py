import numpy as np
from string import punctuation
from copy import copy
from collections import defaultdict
import operator as op
from bisect import bisect
from pyteomics import parser, mass, fasta, auxiliary as aux, mgf, mzml
from . import scoring, utils
try:
    from pyteomics import cmass
except ImportError:
    cmass = mass

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
    for spec in utils.iterate_spectra(fname):
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

    fcharge = settings.getint('scoring', 'maximum fragment charge')
    return {'rel': rel, 'aa_mass': aa_mass, 'acc_l': acc_l, 'acc_r': acc_r, 'acc_frag': acc_frag,
            'unit': unit, 'nmods': nmods, 'maxmods': maxmods, 'fragcharge': fcharge, 'settings': settings}

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
    shifts_and_pime = utils.get_shifts_and_pime(settings)
    dm_l = acc_l * m / 1.0e6 if rel else acc_l# * c FIXME
    dm_r = acc_r * m / 1.0e6 if rel else acc_r# * c FIXME
    cand_idx = []
    for i in shifts_and_pime:
        start = nmasses.searchsorted(m + i - dm_l)
        end = nmasses.searchsorted(m + i + dm_r)       
        if start != end:
            cand_idx.extend(idx[start:end])
    cand_idx = list(set(cand_idx))
    cand_spectra = spectra[cand_idx]
    if settings.has_option('scoring', 'condition'):
        cond = settings.get('scoring', 'condition')
    else:
        cond = None
    if isinstance(cond, str) and cond.strip():
        cond = utils.import_(cond)
    if cond:
        cand_spectra = [c for c in cand_spectra if cond(c, seqm, settings)]
    theor = utils.theor_spectrum(seqm, maxcharge=kwargs['fragcharge'], aa_mass=kwargs['aa_mass'])
    results = [scoring._hyperscore(copy(s), theor, kwargs['acc_frag']) for s in cand_spectra] # FIXME (use score from settings?)
    
    results = [(x.pop('score'), utils.get_title(s), x, m) for x, s in zip(results, cand_spectra)]
    results.sort(reverse=True)
    # results = np.array(results, dtype=[('score', np.float32), ('title', np.str_, 30), ('spectrum', np.object_), ('info', np.object_)])
    return peptide, results

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
            c[-1][4]['mzdiff'] = {'Da': s['nm'][nidx] - pnm}
            c[-1][4]['mzdiff']['ppm'] = 1e6 * c[-1][4]['mzdiff']['Da'] / pnm
            evalues.append(-1./score if -score else 1e6)
        c = np.array(c, dtype=dtype)
        yield {'spectrum': s, 'candidates': c, 'e-values': evalues}

    # return spec_results
    # return (func(p) for p in peps)
    # return utils.multimap(n, func, peps)


