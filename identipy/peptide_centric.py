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
from copy import copy
try:
    from pyteomics import cmass
except ImportError:
#   logger.warning('cmass could not be imported')
    cmass = mass

# try:
# import pyximport; pyximport.install()
from .cutils import theor_spectrum
#from .utils import theor_spectrum
# except:
#     logger.info('Cython modules were not loaded...')
#     from .utils import theor_spectrum
from .utils import reshape_theor_spectrum
# from .scoring import RNHS_ultrafast
from .cutils import RNHS_ultrafast

def prepare_peptide_processor(fname, settings):

    global spectra
    global nmasses
    global nmasses_set
    global titles
    global t2s
    global charges
    global effcharges
    # global best_res
    global fulls_global
    # spectra = {}
    # titles = {}
    # best_res = {}
    # nmasses = {}
    # t2s = {}
    # charges = {}
    spectra = []
    titles = []
    nmasses = []
    charges = []
    effcharges = []
    fulls_global = {}

    nmasses_set = set()

    try:
        fast_first_stage = settings.getint('misc', 'fast first stage')
    except:
        fast_first_stage = 0

    t2s = {}
    best_res = {}
    maxcharges = {}
    fcharge = settings.getint('scoring', 'maximum fragment charge')
    ch_range = range(settings.getint('search', 'minimum charge'),
                1 + settings.getint('search', 'maximum charge'))
    # if fast_first_stage:
    #     fcharge = 1
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
    if not ptol_unit == 'Da' or prec_acc_Da < 1.0:
        prec_acc_Da = False
    else:
        prec_acc_Da = prec_acc_Da
    #     prec_acc_Da = prec_acc_Da * 1e-6 * 5000

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
                    # nmasses.setdefault(effc, []).append(m)
                    # spectra.setdefault(effc, []).append(ps)
                    # titles.setdefault(effc, []).append(ttl)
                    # charges.setdefault(effc, []).append(c)
                    nmasses.append(m)
                    spectra.append(ps)
                    titles.append(ttl)
                    charges.append(c)
                    effcharges.append(effc)
                    ps.setdefault('nm', {})[c] = m
        num_spectra = len(spectra)
        # num_spectra = sum(map(len, spectra.itervalues()))
        logger.info('%s spectra pass quality criteria.', num_spectra)

            # for specs, ttls in zip(spectra[c], titles[c]):
            #     for tmpval in specs['idict']:
            #         fulls_global[tmpval].append(ttls)

        # logger.info(list(fulls_global.items())[:1])

    else:
        num_spectra = len(spectra)
        # num_spectra = sum(map(len, spectra.itervalues()))
        logger.info('Reusing %s spectra from previous run.', num_spectra)

    fulls_global = {}
    i = np.argsort(nmasses)
    nmasses = np.array(nmasses)[i]
    spectra = np.array(spectra)[i]
    titles = np.array(titles)[i]
    charges = np.array(charges)[i]
    effcharges = np.array(effcharges)[i]

    if not ptol_unit == 'Da':
        max_prec_acc_Da = max(nmasses) * 1e-6 * max(abs(lptol), abs(rptol))
    else:
        max_prec_acc_Da = max(abs(lptol), abs(rptol))
    tmp = (nmasses / max_prec_acc_Da).astype(int)
    nmasses_set.update(tmp)
    nmasses_set.update(tmp+1)
    nmasses_set.update(tmp-1)

    if prec_acc_Da:
        nmasses_conv = nmasses / prec_acc_Da
        nmasses_conv = nmasses_conv.astype(int)

        tmp_dict = {}
        for idx, nm in enumerate(nmasses_conv):
            if nm not in tmp_dict:
                tmp_dict[nm] = {}
            if nm+1 not in tmp_dict:
                tmp_dict[nm+1] = {}
            if nm-1 not in tmp_dict:
                tmp_dict[nm-1] = {}
            for spval in spectra[idx]['idict']:
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

        del nmasses_conv

        fulls_global = tmp_dict

    # fulls_global = {}
    # for c in list(spectra):
    #     i = np.argsort(nmasses[c])
    #     nmasses[c] = np.array(nmasses[c])[i]
    #     spectra[c] = np.array(spectra[c])[i]
    #     titles[c] = np.array(titles[c])[i]
    #     charges[c] = np.array(charges[c])[i]

    #     if prec_acc_Da:
    #         nmasses_conv = nmasses[c] / prec_acc_Da
    #         nmasses_conv = nmasses_conv.astype(int)

    #         tmp_dict = {}
    #         for idx, nm in enumerate(nmasses_conv):
    #             if nm not in tmp_dict:
    #                 tmp_dict[nm] = {}#defaultdict(list)
    #             if nm+1 not in tmp_dict:
    #                 tmp_dict[nm+1] = {}#defaultdict(list)
    #             if nm-1 not in tmp_dict:
    #                 tmp_dict[nm-1] = {}#defaultdict(list)
    #             for spval in spectra[c][idx]['idict']:
    #                 # tmp_dict[nm][spval].append(idx)
    #                 # tmp_dict[nm+1][spval].append(idx)
    #                 # tmp_dict[nm-1][spval].append(idx)
    #                 if spval not in tmp_dict[nm]:
    #                     tmp_dict[nm][spval] = [idx, ]
    #                 else:
    #                     tmp_dict[nm][spval].append(idx)
    #                 if spval not in tmp_dict[nm+1]:
    #                     tmp_dict[nm+1][spval] = [idx, ]
    #                 else:
    #                     tmp_dict[nm+1][spval].append(idx)
    #                 if spval not in tmp_dict[nm-1]:
    #                     tmp_dict[nm-1][spval] = [idx, ]
    #                 else:
    #                     tmp_dict[nm-1][spval].append(idx)

    #         del nmasses_conv

    #         fulls_global[c] = tmp_dict

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
                from cutils import RNHS_fast_basic as score_fast_basic
            except: 
                score_fast = utils.import_(settings.get('scoring', 'score') + '_fast')
                score_fast = utils.import_(settings.get('scoring', 'score') + '_fast_basic')
        else:
            score_fast = utils.import_(settings.get('scoring', 'score') + '_fast')
            score_fast = utils.import_(settings.get('scoring', 'score') + '_fast_basic')
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
            'cond': cond, 'score': score, 'score_fast': score_fast, 'score_fast_basic': score_fast_basic,
            'settings': settings, 'max_v': num_spectra, 'prec_acc_Da': prec_acc_Da, 'max_prec_acc_Da': max_prec_acc_Da}

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

    max_prec_acc_Da = kwargs.get('max_prec_acc_Da')

    rel = kwargs['rel']
    acc_l = kwargs['acc_l']
    acc_r = kwargs['acc_r']
    settings = kwargs['settings']

    use_allowed_ions = kwargs.get('use_allowed_ions')
    allowed_ions = kwargs.get('allowed_ions')
    if kwargs['fast first stage']:
        use_allowed_ions = 1
        allowed_ions = {('b', 1), ('y', 1), ('a', 1), ('c', 1), ('x', 1), ('y-17', 1)}
    bions_map = kwargs.get('rank_map')
    yions_map = kwargs.get('rank_map_unf')
    # xions_map = kwargs.get('rank_map_unf_new')

    # logger.info(bions_map)
    shifts_and_pime = kwargs['sapime']
    theor = {}
    theoretical_set = {}
    cand_idx = {}
    stored_value = False
    if rel:
        dm_l = acc_l * m / 1.0e6
        dm_r = acc_r * m / 1.0e6
    elif not rel:
        dm_l = acc_l
        dm_r = acc_r
    # for c in spectra:

    idx = set()
    for shift in shifts_and_pime:
        if int((m + shift)/max_prec_acc_Da) in nmasses_set:
            start = nmasses.searchsorted(m + shift - dm_l)
            end = nmasses.searchsorted(m + shift + dm_r)
            if end - start:
                idx.update(range(start, end))
    if kwargs['cond']:
        idx2 = set()
        for i in idx:
            cond_val, stored_value = kwargs['cond'](spectra[i], seqm, settings, stored_value)
            if cond_val:
                idx2.add(i)
        idx = idx2

    if idx:
        cand_idx = idx
        reshaped = {}
        for c in set(effcharges[i] for i in idx):
            # theor[c], theoretical_set[c] = theor_spectrum(seqm, maxcharge=c, aa_mass=kwargs['aa_mass'], reshape=False,
            #                                                 acc_frag=kwargs['acc_frag'], nterm_mass = nterm_mass,
            #                                                 cterm_mass = cterm_mass, nm=m, simple=kwargs['fast first stage'])
            theor[c], theoretical_set[c] = theor_spectrum(seqm, maxcharge=c, aa_mass=kwargs['aa_mass'], reshape=False,
                                                            acc_frag=kwargs['acc_frag'], nterm_mass = nterm_mass,
                                                            cterm_mass = cterm_mass, nm=m, allowed_ions=allowed_ions, use_allowed_ions=use_allowed_ions)
            reshaped[c] = False
        # reshaped = False

    results = []
    # for ind in cand_idx:
    ind = cand_idx
    # reshaped = False
    if kwargs['prec_acc_Da']:
        print('HERERERE')
        fulls_global_charge = fulls_global
        nm_key = int(m / kwargs['prec_acc_Da'])
        cur_idict = fulls_global_charge.get(nm_key, dict())
        fc_max = max(theor.keys())
        idx_new = RNHS_ultrafast(cur_idict, theoretical_set[fc_max], kwargs['min_matched'], best_res, ind, kwargs['max_v'])
    else:
        idx_new = ind
    if idx_new:
        # logger.info(len(idx_new))
        for i in idx_new:
            # st = utils.get_title(s)
            # if idx_new.count(st) >= kwargs['min_matched']:#st in idx_new:
            # if i in idx_new:
            fc = effcharges[i]
            s = spectra[i]
            st = titles[i]
            chim = ('params' in s and 'isowidthdiff' in s['params'] and abs(float(s['params']['isowidthdiff'])) >= 0.1)
            spcharge = charges[i]
            # neutral_mass, charge_state, RT = get_info(res['spectrum'], res, settings, aa_mass)
            if kwargs['score_fast']:
                if 1:#-idx_new[st] <= best_res.get(st, 0):
                    # logger.info((-idx_new[st], best_res.get(st, 0)))
                    if not bions_map:
                        # print('HERE1')
                        hf = kwargs['score_fast_basic'](s['fastset'], s['idict'], theoretical_set[fc], kwargs['min_matched'])
                    else:
                        # print('HERE2')
                        # hf = kwargs['score_fast_basic'](s['fastset'], s['idict'], theoretical_set[fc], kwargs['min_matched'])
                        hf = kwargs['score_fast'](s['fastset'], s['idict'], theoretical_set[fc], kwargs['min_matched'], bions_map[chim][spcharge])
                    if hf[0]:
                        if -hf[1] <= best_res.get(st, 0):
                            if kwargs['fast first stage']:
                                sc = hf[1]
                                score = {'match': [], 'sumI': 1, 'dist': [], 'total_matched': 999, 'score_std': 0}
                            else:
                                if not reshaped[fc]:
                                    theor[fc] = reshape_theor_spectrum(theor[fc])
                                    reshaped[fc] = True
                                score = kwargs['score'](s, theor[fc], kwargs['acc_frag'], kwargs['acc_frag_ppm'], position=aachange_pos, bions_map=bions_map[chim][spcharge], yions_map=yions_map[chim][spcharge])#, xions_map=xions_map[chim][spcharge])#settings.getfloat('search', 'product accuracy ppm'))  # FIXME (?)
                                sc = score.pop('score')
                            if -sc <= best_res.get(st, 0) and score.pop('total_matched') >= kwargs['min_matched']:
                                results.append((sc, st, charges[i], score))
            else:
                if not reshaped[fc]:
                    theor[fc] = reshape_theor_spectrum(theor[fc])
                    reshaped[fc] = True
                score = kwargs['score'](s, theor[fc], kwargs['acc_frag'], kwargs['acc_frag_ppm'], position=aachange_pos, bions_map=bions_map[chim][spcharge], yions_map=yions_map[chim][spcharge])#, xions_map=xions_map[chim][spcharge])#settings.getfloat('search', 'product accuracy ppm'))  # FIXME (?)
                sc = score.pop('score')
                if -sc <= best_res.get(st, 0) and score.pop('total_matched') >= kwargs['min_matched']:
                    results.append((sc, st, charges[i], score))

    if results:
        return seqm, m, snp_label, results


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

    try:
        kwargs['rank_map'] = settings.get('search', 'rank_map')
        kwargs['rank_map_unf'] = settings.get('search', 'rank_map_unf')
        # kwargs['rank_map_unf_new'] = settings.get('search', 'rank_map_unf_new')
        kwargs['allowed_ions'] = settings.get('search', 'allowed_ions')
        kwargs['use_allowed_ions'] = 1
    except:
        kwargs['rank_map'] = False
        kwargs['rank_map_unf'] = False
        # kwargs['rank_map_unf_new'] = False
        kwargs['allowed_ions'] = set()
        kwargs['use_allowed_ions'] = 0

    logger.info('Running the search ...')
    n = settings.getint('performance', 'processes')
    leg = {}
    if settings.has_option('misc', 'legend'):
        leg = settings.get('misc', 'legend')

    try:
        kwargs['best_peptides'] = settings.get('scoring', 'best peptides')
    except:
        kwargs['best_peptides'] = False

    best_res_raw, best_res = utils.multimap(n, func, peps, **kwargs)

    for spec_t, v in best_res_raw.items():
        peptide, m, snp_label, score, st, c, info = v
        spec_results[spec_t]['spectrum'] = t2s[spec_t]
        info['pep_nm'] = m
        info['charge'] = c
        spec_results[spec_t]['top_scores'] = -score
        spec_results[spec_t]['sequences'] = peptide
        spec_results[spec_t]['info'] = info
        spec_results[spec_t]['snp_label'] = snp_label

    maxlen = settings.getint('search', 'peptide maximum length')
    dtype = np.dtype([('score', np.float64),
        ('seq', np.str_, maxlen), ('note', np.str_, 1),
        ('charge', np.int8), ('info', np.object_), ('sumI', np.float64), ('fragmentMT', np.float64), ('snp_label', np.str_, 15), ('nextscore_std', np.float64)])
    for spec_name, val in spec_results.iteritems():
        s = val['spectrum']
        c = []
        evalues = []
        score = val['top_scores']
        mseq = val['sequences']
        seq = mseq
        info = val['info']
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

