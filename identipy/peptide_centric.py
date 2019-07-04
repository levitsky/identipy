import numpy as np
from string import punctuation
from collections import defaultdict
import operator as op
import random
from bisect import bisect
from pyteomics import parser, mass, fasta, auxiliary as aux, mgf, mzml
from . import scoring, utils
import logging
logger = logging.getLogger(__name__)

try:
    from pyteomics import cmass
except ImportError:
#   logger.warning('cmass could not be imported')
    cmass = mass

# try:
# import pyximport; pyximport.install()
from .cutils import theor_spectrum
# from utils import theor_spectrum
# except:
#     logger.info('Cython modules were not loaded...')
#     from .utils import theor_spectrum
from .utils import reshape_theor_spectrum
# from .scoring import RNHS_ultrafast
from .cutils import RNHS_ultrafast

spectra = {}
titles = {}
best_res = {}
nmasses = {}
t2s = {}
charges = {}
fulls_global = {}
def prepare_peptide_processor(fname, settings):
    global spectra
    global nmasses
    global titles
    global t2s
    global charges
    # global best_res
    global fulls_global
    best_res = {}
    maxcharges = {}
    fcharge = settings.getint('scoring', 'maximum fragment charge')
    ch_range = range(settings.getint('search', 'minimum charge'),
                1 + settings.getint('search', 'maximum charge'))
    for c in ch_range:
        maxcharges[c] = max(1, min(fcharge, c-1) if fcharge else c-1)


    params = {}
    params['maxpeaks'] = settings.getint('scoring', 'maximum peaks')
    params['minpeaks'] = settings.getint('scoring', 'minimum peaks')
    params['dynrange'] = settings.getfloat('scoring', 'dynamic range')
    params['acc'] = settings.getfloat('search', 'product accuracy')
    params['min_mz'] = settings.getfloat('search', 'product minimum m/z')
    params.update(utils._charge_params(settings)) 
    params['dacc'] = settings.getfloat('input', 'deisotoping mass tolerance')
    params['deisotope'] = settings.getboolean('input', 'deisotope')
    params['tags'] = utils.get_tags(settings.get('output', 'tags'))
    rapid_check = settings.getint('search', 'rapid_check')
    
    ptol_unit = settings.get('search', 'precursor accuracy unit')
    lptol = settings.getfloat('search', 'precursor accuracy left')
    rptol = settings.getfloat('search', 'precursor accuracy right')
    prec_acc_Da = max(abs(lptol), abs(rptol))
    # if not ptol_unit == 'Da' or prec_acc_Da < 0.1:
    #     prec_acc_Da = False
    prec_acc_Da = prec_acc_Da
    if not ptol_unit == 'Da':
        prec_acc_Da = prec_acc_Da * 1e-6 * 5000

    if not spectra:
        logger.info('Reading spectra ...')
        if not rapid_check:
            tmp_spec = utils.iterate_spectra(fname)
        else:
            tmp_spec = [spec for spec in utils.iterate_spectra(fname)]
            if len(tmp_spec) >= 2000:
                tmp_spec = random.sample(tmp_spec, 2000)
        for spec in tmp_spec:
            ps = utils.preprocess_spectrum(spec, params)
            if ps is not None:
                ttl = utils.get_title(ps)
                t2s[ttl] = ps
                for m, c in utils.neutral_masses(ps, params):
                    effc = maxcharges[c]
                    nmasses.setdefault(effc, []).append(m)
                    spectra.setdefault(effc, []).append(ps)
                    titles.setdefault(effc, []).append(ttl)
                    charges.setdefault(effc, []).append(c)
                    ps.setdefault('nm', {})[c] = m
        num_spectra = sum(map(len, spectra.itervalues()))
        logger.info('%s spectra pass quality criteria.', num_spectra)

            # for specs, ttls in zip(spectra[c], titles[c]):
            #     for tmpval in specs['idict']:
            #         fulls_global[tmpval].append(ttls)

        # logger.info(list(fulls_global.items())[:1])

    else:
        num_spectra = sum(map(len, spectra.itervalues()))
        logger.info('Reusing %s spectra from previous run.', num_spectra)

    fulls_global = {}
    for c in list(spectra):
        i = np.argsort(nmasses[c])
        nmasses[c] = np.array(nmasses[c])[i]
        spectra[c] = np.array(spectra[c])[i]
        titles[c] = np.array(titles[c])[i]
        charges[c] = np.array(charges[c])[i]

        nmasses_conv = nmasses[c] / prec_acc_Da
        nmasses_conv = nmasses_conv.astype(int)

        tmp_dict = {}
        for idx, nm in enumerate(nmasses_conv):
            if nm not in tmp_dict:
                tmp_dict[nm] = {}#defaultdict(list)
            if nm+1 not in tmp_dict:
                tmp_dict[nm+1] = {}#defaultdict(list)
            if nm-1 not in tmp_dict:
                tmp_dict[nm-1] = {}#defaultdict(list)
            for spval in spectra[c][idx]['idict']:
                # tmp_dict[nm][spval].append(idx)
                # tmp_dict[nm+1][spval].append(idx)
                # tmp_dict[nm-1][spval].append(idx)
                if spval not in tmp_dict[nm]:
                    tmp_dict[nm][spval] = [idx, ]
                else:
                    tmp_dict[nm][spval].append(idx)
                if spval not in tmp_dict[nm+1]:
                    tmp_dict[nm+1][spval] = [idx, ]
                else:
                    tmp_dict[nm+1][spval].append(idx)
                if spval not in tmp_dict[nm-1]:
                    tmp_dict[nm-1][spval] = [idx, ]
                else:
                    tmp_dict[nm-1][spval].append(idx)

        fulls_global[c] = tmp_dict

    utils.set_mod_dict(settings)

    mods = settings.get('modifications', 'variable')
    maxmods = settings.getint('modifications', 'maximum variable mods')
    leg = settings.get('misc', 'legend')
    punct = set(punctuation)
    nmods = [(p, mod[1], mod[2]) for p, mod in leg.iteritems() if p in punct]

    aa_mass = utils.get_aa_mass(settings)
    score = utils.import_(settings.get('scoring', 'score'))
    try:
        score_fast_name = settings.get('scoring', 'score') + '_fast'
        if score_fast_name == 'identipy.scoring.RNHS_fast':
            try:
                from cutils import RNHS_fast as score_fast
            except: 
                score_fast = utils.import_(settings.get('scoring', 'score') + '_fast')
        else:
            score_fast = utils.import_(settings.get('scoring', 'score') + '_fast')
    except Exception as e:
        score_fast = False
        logging.debug('No fast score imported: %s', e)
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
            'settings': settings, 'max_v': num_spectra, 'prec_acc_Da': prec_acc_Da}

def peptide_processor_iter_isoforms(peptide, best_res, **kwargs):
    nmods, maxmods = op.itemgetter('nmods', 'maxmods')(kwargs)
    if nmods and maxmods:
        out = []
        for form in utils.custom_isoforms(peptide, variable_mods=nmods, maxmods=maxmods, snp=kwargs['snp']):
            res = peptide_processor(form, best_res, **kwargs)
            if res:
                out.append(res)
        if out:
            return out
    else:
        res = peptide_processor(peptide, best_res, **kwargs)
        if res:
            return [res, ]


def peptide_processor(peptide, best_res, **kwargs):
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
        aachange_pos = False
    else:
        seqm = peptide
        aachange_pos = False
        snp_label = False

    nterm_mass = kwargs.get('nterm_mass')
    cterm_mass = kwargs.get('cterm_mass')
    m = utils.custom_mass(seqm, aa_mass=kwargs['aa_mass'], nterm_mass = nterm_mass, cterm_mass = cterm_mass)
    # m = cmass.fast_mass(seqm, aa_mass=kwargs['aa_mass']) + (nterm_mass - 1.007825) + (cterm_mass - 17.002735)
    rel = kwargs['rel']
    acc_l = kwargs['acc_l']
    acc_r = kwargs['acc_r']
    settings = kwargs['settings']
    shifts_and_pime = kwargs['sapime']
    theor = {}
    theoretical_set = {}
    cand_idx = {}
    stored_value = False
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
            if end - start:
                idx.update(range(start, end))
        if kwargs['cond']:
            idx2 = set()
            for i in idx:
                cond_val, stored_value = kwargs['cond'](spectra[c][i], seqm, settings, stored_value)
                if cond_val:
                    idx2.add(i)
            idx = idx2

        if idx:
            cand_idx[c] = idx
            theor[c], theoretical_set[c] = theor_spectrum(seqm, maxcharge=c, aa_mass=kwargs['aa_mass'], reshape=False,
                                                          acc_frag=kwargs['acc_frag'], nterm_mass = nterm_mass,
                                                          cterm_mass = cterm_mass, nm=m)
            reshaped = False

    results = []
    for fc, ind in cand_idx.iteritems():
        reshaped = False
        if kwargs['prec_acc_Da']:
            fulls_global_charge = fulls_global[fc]
            nm_key = int(m / kwargs['prec_acc_Da'])
            cur_idict = fulls_global_charge.get(nm_key, dict())
            idx_new = RNHS_ultrafast(cur_idict, theoretical_set[fc], kwargs['min_matched'], best_res, ind, kwargs['max_v'])
        else:
            idx_new = idx
        if idx_new:
        # logger.info(len(idx_new))
            for i in idx_new:
                # st = utils.get_title(s)
                # if idx_new.count(st) >= kwargs['min_matched']:#st in idx_new:
                # if i in idx_new:
                s = spectra[fc][i]
                st = titles[fc][i]
                if kwargs['score_fast']:
                    if 1:#-idx_new[st] <= best_res.get(st, 0):
                        # logger.info((-idx_new[st], best_res.get(st, 0)))
                        hf = kwargs['score_fast'](s['fastset'], s['idict'], theoretical_set[fc], kwargs['min_matched'])
                        if hf[0]:
                            if -hf[1] <= best_res.get(st, 0):
                                if kwargs['fast first stage']:
                                    sc = hf[1]
                                    score = {'match': [], 'sumI': 1, 'dist': [], 'total_matched': 999, 'score_std': 0}
                                else:
                                    if not reshaped:
                                        theor[fc] = reshape_theor_spectrum(theor[fc])
                                        reshaped = True
                                    score = kwargs['score'](s, theor[fc], kwargs['acc_frag'], kwargs['acc_frag_ppm'], position=aachange_pos)#settings.getfloat('search', 'product accuracy ppm'))  # FIXME (?)
                                    sc = score.pop('score')
                                if -sc <= best_res.get(st, 0) and score.pop('total_matched') >= kwargs['min_matched']:
                                    # results.append((sc, st, score, m, charges[fc][i], snp_label))
                                    # results.append((sc, st, score, charges[fc][i]))
                                    # results.append((sc, fc, i))
                                    results.append((sc, st, charges[fc][i], score))
                    # else:
                    #     logger.info('FAIL')
                else:
                    # st = utils.get_title(s)
                    if not reshaped:
                        theor[fc] = reshape_theor_spectrum(theor[fc])
                        reshaped = True
                    score = kwargs['score'](s, theor[fc], kwargs['acc_frag'], kwargs['acc_frag_ppm'], position=aachange_pos)#settings.getfloat('search', 'product accuracy ppm'))  # FIXME (?)
                    sc = score.pop('score')
                    if -sc <= best_res.get(st, 0) and score.pop('total_matched') >= kwargs['min_matched']:
                        # results.append((sc, st, score, m, charges[fc][i], snp_label))
                        # results.append((sc, st, score, charges[fc][i]))
                        # results.append((sc, fc, i))
                        results.append((sc, st, charges[fc][i], score))


    # results.sort(reverse=True, key=op.itemgetter(0))
    # results = np.array(results, dtype=[('score', np.float32), ('title', np.str_, 30), ('spectrum', np.object_), ('info', np.object_)])
    if results:
        return seqm, m, snp_label, results
    # return seqm, []

def recalc_best_candidates(seqm, m, fc, i, kwargs):
    theor, _ = theor_spectrum(seqm, maxcharge=fc, aa_mass=kwargs['aa_mass'], reshape=True,
                                                            acc_frag=kwargs['acc_frag'], nterm_mass = kwargs['nterm_mass'],
                                                            cterm_mass = kwargs['cterm_mass'], nm=m)
    s = spectra[fc][i]
    score = kwargs['score'](s, theor, kwargs['acc_frag'], kwargs['acc_frag_ppm'], position=False)#settings.getfloat('search', 'product accuracy ppm'))  # FIXME (?)
    return score

def process_peptides(fname, settings):
    spec_results = defaultdict(dict)
    peps = utils.peptide_gen(settings)
    kwargs = prepare_peptide_processor(fname, settings)
    func = peptide_processor_iter_isoforms
    kwargs['min_matched'] = settings.getint('output', 'minimum matched')
    kwargs['snp'] = settings.getint('search', 'snp')
    kwargs['nterm_mass'] = settings.getfloat('modifications', 'protein nterm cleavage')
    kwargs['cterm_mass'] = settings.getfloat('modifications', 'protein cterm cleavage')
    kwargs['qsize'] = settings.getint('performance', 'out queue size')
    logger.info('Running the search ...')
    n = settings.getint('performance', 'processes')
    leg = {}
    if settings.has_option('misc', 'legend'):
        leg = settings.get('misc', 'legend')




    # rel = kwargs['rel']
    # acc_l = kwargs['acc_l']
    # acc_r = kwargs['acc_r']
    # if not rel and abs(acc_r + acc_l) >= 10:
    #     logger.info('Run closed search')
    #     kwargs['acc_l'] = 20
    #     kwargs['acc_r'] = 20
    #     kwargs['rel'] = True


    #     best_res_raw, best_res = utils.multimap(n, func, peps, **kwargs)
    #     # for y in utils.multimap(n, func, peps, **kwargs):
    #     #     for x in y:
    #     #         # if x[3] is not None:
    #     #         peptide, m, snp_label, result = x

    #     #         for score, fc, i, info in result:
    #     #             spec_t = titles[fc][i]
    #     #             if -score <= best_res.get(spec_t, 0):
    #     #                 best_res_raw[spec_t] = [peptide, m, snp_label, score, fc, i, info]
    #     #                 best_res[spec_t] = -score

    #     peps = utils.peptide_gen(settings, clear_seen_peptides=True)
    #     kwargs['acc_l'] = acc_l
    #     kwargs['acc_r'] = acc_r
    #     kwargs['rel'] = rel
    #     logger.info('Closed search is finished')




    # for y in utils.multimap(n, func, peps, **kwargs):
    #     for x in y:
    #         # if x[3] is not None:
    #         peptide, m, snp_label, result = x

    #         for score, fc, i, info in result:
    #             spec_t = titles[fc][i]
    #             if -score <= best_res.get(spec_t, 0):
    #                 best_res_raw[spec_t] = [peptide, m, snp_label, score, fc, i, info]
    #                 best_res[spec_t] = -score    

    best_res_raw, best_res = utils.multimap(n, func, peps, **kwargs)
    # best_res_raw, best_res = utils.multimap(n, func, peps, best_res, best_res_raw, **kwargs)
    # for y in utils.multimap(n, func, peps, best_res, best_res_raw, **kwargs):
    #     for x in y:
    #         # if x[3] is not None:
    #         peptide, m, snp_label, result = x

    #         for score, fc, i, info in result:
    #             spec_t = titles[fc][i]
    #             if -score <= best_res.get(spec_t, 0):
    #                 best_res_raw[spec_t] = [peptide, m, snp_label, score, fc, i, info]
    #                 best_res[spec_t] = -score    

    for spec_t, v in best_res_raw.items():
        peptide, m, snp_label, score, st, c, info = v
        spec_results[spec_t]['spectrum'] = t2s[spec_t]
        info['pep_nm'] = m
        info['charge'] = c
        spec_results[spec_t]['top_scores'] = -score
        spec_results[spec_t]['sequences'] = peptide
        spec_results[spec_t]['info'] = info
        spec_results[spec_t]['snp_label'] = snp_label


            # for score, fc, i in result:
            #     spec_t = titles[fc][i]
            #     spec_results[spec_t]['spectrum'] = t2s[spec_t]
            #     top_scores = spec_results[spec_t].setdefault('top_scores', 0)
            #     if -score <= top_scores:
            #         c = charges[fc][i]
            #         best_res[spec_t] = -score
            #         info = recalc_best_candidates(peptide, m, fc, i, kwargs)
            #         info['pep_nm'] = m
            #         info['charge'] = c
            #         spec_results[spec_t]['top_scores'] = -score
            #         spec_results[spec_t]['sequences'] = peptide
            #         spec_results[spec_t]['info'] = info
            #         spec_results[spec_t]['snp_label'] = snp_label

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
        ('charge', np.int8), ('info', np.object_), ('sumI', np.float64), ('fragmentMT', np.float64), ('snp_label', np.str_, 15), ('nextscore_std', np.float64)])
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
            repl = leg[x][1]
            if repl == '-':
                repl = ''
            seq = seq.replace(x, repl)
        pnm = info['pep_nm']
        c.append((-score, mseq, 't' if seq in utils.seen_target else 'd',
            info['charge'], info, info.pop('sumI'), np.median(info.pop('dist')), val['snp_label'], info.pop('score_std')))
        c[-1][4]['mzdiff'] = {'Da': s['nm'][info['charge']] - pnm}
        c[-1][4]['mzdiff']['ppm'] = 1e6 * c[-1][4]['mzdiff']['Da'] / pnm
        evalues.append(-1./score if -score else 1e6)
        c = np.array(c, dtype=dtype)
        yield {'spectrum': s, 'candidates': c, 'e-values': evalues}

