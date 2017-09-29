import numpy as np
from string import punctuation
from collections import defaultdict
import operator as op
from bisect import bisect
from pyteomics import parser, mass, fasta, auxiliary as aux, mgf, mzml
from . import scoring, utils
try:
    from pyteomics import cmass
except ImportError:
    print '[ Warning ] cmass could not be imported'
    cmass = mass

try:
    import pyximport; pyximport.install()
    from cutils import theor_spectrum
except:
    from utils import theor_spectrum

def prepare_peptide_processor(fname, settings):
    global spectra
    global nmasses
    global t2s
    global charges
    global best_res
    spectra = {}
    best_res = {}
    nmasses = {}
    t2s = {}
    charges = {}

    maxcharges = {}
    fcharge = settings.getint('scoring', 'maximum fragment charge')
    ch_range = range(settings.getint('search', 'minimum charge'),
                1 + settings.getint('search', 'maximum charge'))
    for c in ch_range:
        maxcharges[c] = max(1, min(fcharge, c-1) if fcharge else c-1)

    print 'Reading spectra ...'

    params = {}
    params['maxpeaks'] = settings.getint('scoring', 'maximum peaks')
    params['minpeaks'] = settings.getint('scoring', 'minimum peaks')
    params['dynrange'] = settings.getfloat('scoring', 'dynamic range')
    params['acc'] = settings.getfloat('search', 'product accuracy')
    params['min_mz'] = settings.getfloat('search', 'product minimum m/z')
    params['maxcharge'] = settings.getint('search', 'maximum charge') or None
    params['mincharge'] = settings.getint('search', 'minimum charge') or None
    if settings.has_option('search', 'minimum unknown charge'):
        params['min_ucharge'] = max(settings.getint('search', 'minimum unknown charge'), params['mincharge'])
    else:
        params['min_ucharge'] = params['mincharge']
    if settings.has_option('search', 'maximum unknown charge'):
        params['max_ucharge'] = min(settings.getint('search', 'maximum unknown charge'), params['maxcharge'])
    else:
        params['max_ucharge'] = params['maxcharge']

    params['dacc'] = settings.getfloat('input', 'deisotoping mass tolerance')
    params['deisotope'] = settings.getboolean('input', 'deisotope')

    for spec in utils.iterate_spectra(fname):
        ps = utils.preprocess_spectrum(spec, params)
        if ps is not None:
            t2s[utils.get_title(ps)] = ps
            for m, c in utils.neutral_masses(ps, params):
                effc = maxcharges[c]
                nmasses.setdefault(effc, []).append(m)
                spectra.setdefault(effc, []).append(ps)
                charges.setdefault(effc, []).append(c)
                ps.setdefault('nm', {})[c] = m
    print sum(map(len, spectra.itervalues())), 'spectra pass quality criteria.'
    for c in list(spectra):
        i = np.argsort(nmasses[c])
        nmasses[c] = np.array(nmasses[c])[i]
        spectra[c] = np.array(spectra[c])[i]
        charges[c] = np.array(charges[c])[i]

    utils.set_mod_dict(settings)

    mods = settings.get('modifications', 'variable')
    maxmods = settings.getint('modifications', 'maximum variable mods')
    leg = settings.get('misc', 'legend')
    punct = set(punctuation)
    nmods = [(p, mod[1], mod[2]) for p, mod in leg.iteritems() if p in punct]

    aa_mass = utils.get_aa_mass(settings)
    score = utils.import_(settings.get('scoring', 'score'))
    try:
        score_fast = utils.import_(settings.get('scoring', 'score') + '_fast')
    except:
        score_fast = False
    acc_l = settings.getfloat('search', 'precursor accuracy left')
    acc_r = settings.getfloat('search', 'precursor accuracy right')
    acc_frag = settings.getfloat('search', 'product accuracy')
    frag_unit = settings.get('search', 'product accuracy unit')
    if frag_unit == 'ppm':
        acc_frag_ppm = settings.getfloat('search', 'product accuracy ppm')
    else:
        acc_frag_ppm = False
    try:
        fast_first_stage = settings.getint('misc', 'fast first stage')
    except:
        fast_first_stage = 0
    unit = settings.get('search', 'precursor accuracy unit')
    rel = utils.relative(unit)

    if settings.has_option('scoring', 'condition'):
        cond = settings.get('scoring', 'condition')
    else:
        cond = None
    if isinstance(cond, str) and cond.strip():
        cond = utils.import_(cond)

    score = utils.import_(settings.get('scoring', 'score'))

    return {'rel': rel, 'aa_mass': aa_mass, 'acc_l': acc_l, 'acc_r': acc_r, 'acc_frag': acc_frag, 'acc_frag_ppm': acc_frag_ppm,
            'unit': unit, 'nmods': nmods, 'maxmods': maxmods, 'fast first stage': fast_first_stage,
            'sapime': utils.get_shifts_and_pime(settings),
            'cond': cond, 'score': score, 'score_fast': score_fast,
            'settings': settings}

def peptide_processor_iter_isoforms(peptide, **kwargs):
    nmods, maxmods = op.itemgetter('nmods', 'maxmods')(kwargs)
    if nmods and maxmods:
        out = []
        for form in utils.custom_isoforms(peptide, variable_mods=nmods, maxmods=maxmods):
            res = peptide_processor(form, **kwargs)
            if res:
                out.append(res)
        if out:
            return out
    else:
        res = peptide_processor(peptide, **kwargs)
        if res:
            return [res, ]


def peptide_processor(peptide, **kwargs):
    if kwargs['snp']:
        if 'snp' not in peptide:
            seqm = peptide
            aachange_pos = False
            snp_label = 'wild'
        else:
            tmp = peptide.split('snp')
            seqm = tmp[0] + tmp[1].split('at')[0].split('to')[-1] + tmp[2]
            aachange_pos = len(tmp[0]) + 1
            snp_label = tmp[1]
    else:
        seqm = peptide
        aachange_pos = False
        snp_label = False
    m = cmass.fast_mass(seqm, aa_mass=kwargs['aa_mass'])
    rel = kwargs['rel']
    acc_l = kwargs['acc_l']
    acc_r = kwargs['acc_r']
    settings = kwargs['settings']
    shifts_and_pime = kwargs['sapime']
    theor = {}
    theoretical_set = {}
    cand_idx = {}
    if rel:
        dm_l = acc_l * m / 1.0e6
        dm_r = acc_r * m / 1.0e6
    for c in spectra:
        if not rel:
            dm_l = acc_l * c
            dm_r = acc_r * c
        idx = set()
        for shift in shifts_and_pime:
            start = nmasses[c].searchsorted(m + shift - dm_l)
            end   = nmasses[c].searchsorted(m + shift + dm_r)
            idx.update(range(start, end))
        if kwargs['cond']:
            idx = {i for i in idx if kwargs['cond'](spectra[c][i], seqm, settings)}

        if idx:
            cand_idx[c] = idx
            theor[c], theoretical_set[c] = theor_spectrum(seqm, maxcharge=c, aa_mass=kwargs['aa_mass'], reshape=True, acc_frag=kwargs['acc_frag'])

    results = []
    for fc, ind in cand_idx.iteritems():
        # breaker = False
        # if len(ind) > 1:
        #     megaspectra = {'m/z array': np.concatenate([spectra[fc][i]['m/z array'] for i in ind]), 'intensity array': np.concatenate([spectra[fc][i]['intensity array'] for i in ind])}
        #     score = kwargs['score'](megaspectra, theor[fc], kwargs['acc_frag'])  # FIXME (?)
        #     sc = score.pop('score')
        #     if score.pop('total_matched') < kwargs['min_matched'] or all(-sc >= best_res.get(utils.get_title(spectra[fc][i]), 0) for i in ind):
        #         breaker = True
        # if not breaker:
        #     for i in ind:
        #         s = spectra[fc][i]
        #         score = kwargs['score'](s, theor[fc], kwargs['acc_frag']) # FIXME (?)
        #         sc = score.pop('score')
        #         st = utils.get_title(s)
        #         if -sc <= best_res.get(st, 0) and score.pop('total_matched') >= kwargs['min_matched']:
        #             results.append((sc, st, score, m, charges[fc][i]))

        for i in ind:
            s = spectra[fc][i]
            if kwargs['score_fast']:
                hf = kwargs['score_fast'](s['fastset'], theoretical_set[fc], kwargs['min_matched'])
                if hf[0]:
                    st = utils.get_title(s)
                    if -hf[1] <= best_res.get(st, 0):
                        if kwargs['fast first stage']:
                            sc = hf[1]
                            score = {'match': [], 'sumI': 1, 'dist': [], 'total_matched': 999}
                        else:
                            score = kwargs['score'](s, theor[fc], kwargs['acc_frag'], kwargs['acc_frag_ppm'], position=aachange_pos)#settings.getfloat('search', 'product accuracy ppm'))  # FIXME (?)
                            sc = score.pop('score')
                        # st = utils.get_title(s)
                        if -sc <= best_res.get(st, 0) and score.pop('total_matched') >= kwargs['min_matched']:
                            results.append((sc, st, score, m, charges[fc][i], snp_label))
            else:
                st = utils.get_title(s)
                score = kwargs['score'](s, theor[fc], kwargs['acc_frag'], kwargs['acc_frag_ppm'], position=aachange_pos)#settings.getfloat('search', 'product accuracy ppm'))  # FIXME (?)
                sc = score.pop('score')
                if -sc <= best_res.get(st, 0) and score.pop('total_matched') >= kwargs['min_matched']:
                    results.append((sc, st, score, m, charges[fc][i], snp_label))


    results.sort(reverse=True, key=op.itemgetter(0))
    # results = np.array(results, dtype=[('score', np.float32), ('title', np.str_, 30), ('spectrum', np.object_), ('info', np.object_)])
    if results:
        return seqm, results

def process_peptides(fname, settings):

    spec_results = defaultdict(dict)
    peps = utils.peptide_gen(settings)
    kwargs = prepare_peptide_processor(fname, settings)
    func = peptide_processor_iter_isoforms
    kwargs['min_matched'] = settings.getint('output', 'minimum matched')
    kwargs['snp'] = settings.getint('search', 'snp')
    print 'Running the search ...'
    n = settings.getint('performance', 'processes')
    leg = {}
    if settings.has_option('misc', 'legend'):
        leg = settings.get('misc', 'legend')
    for y in utils.multimap(n, func, peps, **kwargs):
        for x in y:
            if x is not None:
                peptide, result = x
                for score, spec_t, info, m, c, snp_label in result:
                    spec_results[spec_t]['spectrum'] = t2s[spec_t]
                    top_scores = spec_results[spec_t].setdefault('top_scores', 0)
                    if -score <= top_scores:
                        best_res[spec_t] = -score
                        info['pep_nm'] = m
                        info['charge'] = c
                        spec_results[spec_t]['top_scores'] = -score
                        spec_results[spec_t]['sequences'] = peptide
                        spec_results[spec_t]['info'] = info
                        spec_results[spec_t]['snp_label'] = snp_label

#               spec_results[spec_t].setdefault('scores', []).append(score) FIXME write histogram
#
#                     top_seqs   = spec_results[spec_t].setdefault('sequences', '')
#                     top_info   = spec_results[spec_t].setdefault('info', [])
#
#                     i = bisect(top_scores, -score)
#                     if nc is None or i < nc:
#                         top_scores.insert(i, -score)
#                         top_seqs.insert(i, peptide)
#                         top_info.insert(i, info)
#                         if nc is not None and len(top_scores) > nc:
#                             top_scores.pop()
#                             top_seqs.pop()
#                             top_info.pop()
    maxlen = settings.getint('search', 'peptide maximum length')
    dtype = np.dtype([('score', np.float64),
        ('seq', np.str_, maxlen), ('note', np.str_, 1),
        ('charge', np.int8), ('info', np.object_), ('sumI', np.float64), ('fragmentMT', np.float64), ('snp_label', np.str_, 15)])

    for spec_name, val in spec_results.iteritems():
        s = val['spectrum']
        c = []
        evalues = []
        score = val['top_scores']
        # for idx, score in enumerate(val['top_scores']):
        mseq = val['sequences']#[idx]
        seq = mseq
        info = val['info']#[idx]
        for x in set(mseq).intersection(punctuation):
            seq = seq.replace(x, leg[x][1])
        pnm = info['pep_nm']
        c.append((-score, mseq, 't' if seq in utils.seen_target else 'd',
            info['charge'], info, info.pop('sumI'), np.median(info.pop('dist')), val['snp_label']))
        c[-1][4]['mzdiff'] = {'Da': s['nm'][info['charge']] - pnm}
        c[-1][4]['mzdiff']['ppm'] = 1e6 * c[-1][4]['mzdiff']['Da'] / pnm
        evalues.append(-1./score if -score else 1e6)
        c = np.array(c, dtype=dtype)
        yield {'spectrum': s, 'candidates': c, 'e-values': evalues}

