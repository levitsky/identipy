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
    spectra = {}
    nmasses = {}
    print 'Reading spectra ...'
    for spec in utils.iterate_spectra(fname):
        ps = utils.preprocess_spectrum(spec, settings)
        if ps is not None:
            for m, c in utils.neutral_masses(ps, settings):
                nmasses.setdefault(c, []).append(m)
                spectra.setdefault(c, []).append(ps)
                ps.setdefault('nm', []).append(m)
                ps.setdefault('ch', []).append(c)
    print sum(map(len, spectra.itervalues())), 'spectra pass quality criteria.'
    for c in list(spectra):
        i = np.argsort(nmasses[c])
        nmasses[c] = np.array(nmasses[c])[i]
        spectra[c] = np.array(spectra[c])[i]

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
    if settings.has_option('scoring', 'condition'):
        cond = settings.get('scoring', 'condition')
    else:
        cond = None
    if isinstance(cond, str) and cond.strip():
        cond = utils.import_(cond)

    return {'rel': rel, 'aa_mass': aa_mass, 'acc_l': acc_l, 'acc_r': acc_r, 'acc_frag': acc_frag,
            'unit': unit, 'nmods': nmods, 'maxmods': maxmods, 'fragcharge': fcharge,
            'sapime': utils.get_shifts_and_pime(settings),
            'ch_range': range(settings.getint('search', 'minimum charge'),
                1 + settings.getint('search', 'maximum charge')),
            'cond': cond,
            'settings': settings}

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
    shifts_and_pime = kwargs['sapime']
    theor = {}
    cand_spectra = {}
    for c in kwargs['ch_range']:
        dm_l = acc_l * m / 1.0e6 if rel else acc_l * c
        dm_r = acc_r * m / 1.0e6 if rel else acc_r * c
        idx = set()
        for shift in shifts_and_pime:
            start = nmasses[c].searchsorted(m + shift - dm_l)
            end   = nmasses[c].searchsorted(m + shift + dm_r)
            idx.update(range(start, end))
        if kwargs['cond']:
            idx = {i for i in idx if cond(spectra[c][i], seqm, settings)}

        if idx:
            cand_spectra[c] = spectra[c][list(idx)]
            maxcharge = max(1, min(kwargs['fragcharge'], c-1) if kwargs['fragcharge'] else c-1)
            theor[c] = utils.theor_spectrum(seqm, maxcharge=maxcharge, aa_mass=kwargs['aa_mass'])

    results = []
    for c, cand in cand_spectra.iteritems():
        for s in cand:
            score = scoring._hyperscore(copy(s), theor[c], kwargs['acc_frag']) # FIXME (use score from settings?)
            results.append((score.pop('score'), utils.get_title(s), score, m, s))
    results.sort(reverse=True)
    # results = np.array(results, dtype=[('score', np.float32), ('title', np.str_, 30), ('spectrum', np.object_), ('info', np.object_)])
    return peptide, results

def process_peptides(fname, settings):

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
                for score, spec_t, info, m, spec in result:
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


