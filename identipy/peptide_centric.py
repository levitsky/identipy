import numpy as np
from string import punctuation
from collections import defaultdict
import random
from pyteomics import mass
from . import utils
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
#from .utils import theor_spectrum
# except:
#     logger.info('Cython modules were not loaded...')
#     from .utils import theor_spectrum
from .utils import reshape_theor_spectrum
# from .scoring import RNHS_ultrafast
# from .cutils import RNHS_ultrafast

def prepare_peptide_processor(fname, settings):

    global_data = list()
    n_proc = utils.get_nprocesses(settings)

    for _ in range(n_proc):
        global_data.append({
            'spectra': [],
            'titles': [],
            'nmasses': [],
            'nmasses_set': set(),
            't2s': {},
            'charges': [],
            'effcharges': [],
            'fulls_global': {},
        })

    logger.debug('global data: %s', len(global_data))

    try:
        fast_first_stage = settings.getint('misc', 'fast first stage')
    except:
        fast_first_stage = 0

    # t2s = {}
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
    params['maxcharges'] = maxcharges
    rapid_check = settings.getint('search', 'rapid_check')

    ptol_unit = settings.get('search', 'precursor accuracy unit')
    lptol = settings.getfloat('search', 'precursor accuracy left')
    rptol = settings.getfloat('search', 'precursor accuracy right')
    prec_acc_Da = False
    # prec_acc_Da = max(abs(lptol), abs(rptol))
    # if ptol_unit != 'Da' or prec_acc_Da < 1.0:
    #     prec_acc_Da = False

    logger.info('Reading spectra ...')
    if not rapid_check:
        tmp_spec = utils.iterate_spectra(fname)
    else:
        tmp_spec = [spec for spec in utils.iterate_spectra(fname)]
        if len(tmp_spec) >= 2000:
            tmp_spec = random.sample(tmp_spec, 2000)

    num_spectra = 0

    tmp_spec2 = []
    nmasses_tmp = []
    charges_tmp = []
    global_data_index_map = {}

    for spec in tmp_spec:
        ps = utils.preprocess_spectrum(spec, params)
        if ps is not None:

            tmp_spec2.append(ps)
            for m, c in utils.neutral_masses(ps, params):
                nmasses_tmp.append(m)
                charges_tmp.append(c)

    nmasses_tmp = np.array(nmasses_tmp)
    idx_t = np.argsort(nmasses_tmp)
    max_nmass = nmasses_tmp[idx_t[-1]]
    max_l = int(len(nmasses_tmp)/n_proc)+1
    for idx, k in enumerate(idx_t):
        global_data_index_map[k] = idx // max_l
    logger.debug('nproc: %d, nmasses: %d, max_l: %d, maximum index: %d',
        n_proc, nmasses_tmp.size, max_l, max(global_data_index_map.values()))

    for ps in tmp_spec2:
        # global_data_index = num_spectra % n_proc
        ttl = utils.get_title(ps)
        # t2s[ttl] = ps
        for m, c in utils.neutral_masses(ps, params):
            global_data_index = global_data_index_map[num_spectra]
            effc = maxcharges[c]
            ps.setdefault('nm', {})[c] = m

            global_data[global_data_index]['t2s'][ttl] = ps
            global_data[global_data_index]['nmasses'].append(m)
            global_data[global_data_index]['spectra'].append(ps)
            global_data[global_data_index]['titles'].append(ttl)
            global_data[global_data_index]['charges'].append(c)
            global_data[global_data_index]['effcharges'].append(effc)

            num_spectra += 1
    logger.info('%s spectra pass quality criteria.', num_spectra)

    if ptol_unit != 'Da':
        max_prec_acc_Da = max_nmass * 1e-6 * max(abs(lptol), abs(rptol))
    else:
        max_prec_acc_Da = max(abs(lptol), abs(rptol))


    for global_data_index in range(n_proc):

        i = np.argsort(global_data[global_data_index]['nmasses'])
        global_data[global_data_index]['nmasses'] = np.array(global_data[global_data_index]['nmasses'])[i]
        global_data[global_data_index]['spectra'] = np.array(global_data[global_data_index]['spectra'])[i]
        global_data[global_data_index]['titles'] = np.array(global_data[global_data_index]['titles'])[i]
        global_data[global_data_index]['charges'] = np.array(global_data[global_data_index]['charges'])[i]
        global_data[global_data_index]['effcharges'] = np.array(global_data[global_data_index]['effcharges'])[i]

        tmp = (global_data[global_data_index]['nmasses'] / max_prec_acc_Da).astype(int)
        global_data[global_data_index]['nmasses_set'].update(tmp)
        global_data[global_data_index]['nmasses_set'].update(tmp+1)
        global_data[global_data_index]['nmasses_set'].update(tmp-1)

        if prec_acc_Da:
            nmasses_conv = global_data[global_data_index]['nmasses'] / max_prec_acc_Da
            nmasses_conv = nmasses_conv.astype(int)

            tmp_dict = {}
            for idx, nm in enumerate(nmasses_conv):
                if nm not in tmp_dict:
                    tmp_dict[nm] = {}
                if nm+1 not in tmp_dict:
                    tmp_dict[nm+1] = {}
                if nm-1 not in tmp_dict:
                    tmp_dict[nm-1] = {}
                for spval in global_data[global_data_index]['spectra'][idx]['idict']:
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

            global_data[global_data_index]['nmasses_set'] = tmp_dict

    utils.set_mod_dict(settings)

    aa_mass = utils.get_aa_mass(settings)
    score = utils.import_(settings.get('scoring', 'score'))
    try:
        score_fast_name = settings.get('scoring', 'score') + '_fast'
        logger.debug('Fast score name: %s', score_fast_name)
        if score_fast_name in {'identipy.scoring.RNHS_fast', 'RNHS_fast'}:
            try:
                from .cutils import RNHS_fast as score_fast
                from .cutils import RNHS_fast_basic as score_fast_basic
            except ImportError as e:
                logger.warning('Could not import from cutils: %s', e.args)
                score_fast = utils.import_(settings.get('scoring', 'score') + '_fast')
                score_fast_basic = utils.import_(settings.get('scoring', 'score') + '_fast_basic')
        else:
            score_fast = utils.import_(settings.get('scoring', 'score') + '_fast')
            score_fast_basic = utils.import_(settings.get('scoring', 'score') + '_fast_basic')
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

    return {'rel': rel, 'aa_mass': aa_mass,
            'acc_l': acc_l, 'acc_r': acc_r, 'acc_frag': acc_frag, 'acc_frag_ppm': acc_frag_ppm,
            'unit': unit,  # 'nmods': nmods, 'maxmods': maxmods,
            'fast first stage': fast_first_stage,
            'sapime': utils.get_shifts_and_pime(settings),
            'cond': cond, 'score': score, 'score_fast': score_fast, 'score_fast_basic': score_fast_basic,
            'settings': settings, 'max_v': num_spectra, 'prec_acc_Da': prec_acc_Da, 'max_prec_acc_Da': max_prec_acc_Da}, global_data


def peptide_processor_iter_isoforms(peptide, best_res, global_data_local, **kwargs):
    res = peptide_processor(peptide, best_res, global_data_local, **kwargs)
    if res:
        return [res, ]

    # nmods, maxmods = op.itemgetter('nmods', 'maxmods')(kwargs)
    # if nmods and maxmods:
    #     out = []
    #     for form in utils.custom_isoforms(peptide, variable_mods=nmods, maxmods=maxmods, snp=kwargs['snp']):
    #         res = peptide_processor(form, best_res, global_data_local, **kwargs)
    #         if res:
    #             out.append(res)
    #     if out:
    #         return out
    # else:
    #     res = peptide_processor(peptide, best_res, global_data_local, **kwargs)
    #     if res:
    #         return [res, ]


def peptide_processor(peptide, best_res, global_data_local, **kwargs):
    spectra = global_data_local['spectra']
    titles = global_data_local['titles']
    nmasses = global_data_local['nmasses']
    nmasses_set = global_data_local['nmasses_set']
    t2s = global_data_local['t2s']
    charges = global_data_local['charges']
    effcharges = global_data_local['effcharges']
    fulls_global = global_data_local['fulls_global']
    seqm, aachange_pos, snp_label, m = peptide

    max_prec_acc_Da = kwargs.get('max_prec_acc_Da')

    nterm_mass = kwargs.get('nterm_mass')
    cterm_mass = kwargs.get('cterm_mass')
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
    elif not rel:
        dm_l = acc_l
        dm_r = acc_r
    # for c in spectra:

    idx = set()
    for shift in shifts_and_pime:
        if int((m + shift)/max_prec_acc_Da) in nmasses_set:
            start = nmasses.searchsorted(m + shift - dm_l)
            end = nmasses.searchsorted(m + shift + dm_r, side='right')
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
            theor[c], theoretical_set[c] = theor_spectrum(seqm, maxcharge=c, aa_mass=kwargs['aa_mass'], reshape=False,
                                                            acc_frag=kwargs['acc_frag'], nterm_mass = nterm_mass,
                                                            cterm_mass = cterm_mass, nm=m)
            reshaped[c] = False
        # reshaped = False

    results = []
    # for ind in cand_idx:
    ind = cand_idx
    # reshaped = False
    idx_new = ind
    # if idx_new and kwargs['prec_acc_Da']:
    #     fulls_global_charge = fulls_global
    #     nm_key = int(m / max_prec_acc_Da)
    #     cur_idict = fulls_global_charge.get(nm_key, dict())
    #     fc_max = max(theor.keys())
    #     idx_new = RNHS_ultrafast(cur_idict, theoretical_set[fc_max], kwargs['min_matched'], best_res, ind, kwargs['max_v'])
            
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
                if 1:
                    hf = kwargs['score_fast_basic'](s['fastset'], s['idict'], theoretical_set[fc], kwargs['min_matched'])
                    if hf[0]:
                        if -hf[1] <= best_res.get(st, 0):
                            if kwargs['fast first stage']:
                                sc = hf[1]
                                score = {'match': [], 'sumI': 1, 'dist': [], 'total_matched': 999, 'score_std': 0}
                            else:
                                if not reshaped[fc]:
                                    theor[fc] = reshape_theor_spectrum(theor[fc])
                                    reshaped[fc] = True
                                score = kwargs['score'](s, theor[fc], kwargs['acc_frag'], kwargs['acc_frag_ppm'], position=aachange_pos) # FIXME (?)
                                sc = score.pop('score')
                            if -sc <= best_res.get(st, 0) and score.pop('total_matched') >= kwargs['min_matched']:
                                results.append((sc, st, charges[i], score))
            else:
                if not reshaped[fc]:
                    theor[fc] = reshape_theor_spectrum(theor[fc])
                    reshaped[fc] = True
                score = kwargs['score'](s, theor[fc], kwargs['acc_frag'], kwargs['acc_frag_ppm'], position=aachange_pos) # FIXME (?)
                sc = score.pop('score')
                if -sc <= best_res.get(st, 0) and score.pop('total_matched') >= kwargs['min_matched']:
                    results.append((sc, st, charges[i], score))

    if results:
        return seqm, m, snp_label, results


def process_peptides(fname, settings):
    logger.debug('Started process_peptides.')
    spec_results = defaultdict(dict)
    peps = utils.peptide_isoforms(settings)
    kwargs, global_data = prepare_peptide_processor(fname, settings)
    func = peptide_processor_iter_isoforms
    kwargs['min_matched'] = settings.getint('output', 'minimum matched')
    kwargs['snp'] = settings.getint('search', 'snp')
    kwargs['nterm_mass'] = settings.getfloat('modifications', 'protein nterm cleavage')
    kwargs['cterm_mass'] = settings.getfloat('modifications', 'protein cterm cleavage')
    kwargs['qsize'] = settings.getint('performance', 'out queue size')

    logger.info('Running the search ...')
    n = utils.get_nprocesses(settings)
    leg = {}
    if settings.has_option('misc', 'legend'):
        leg = settings.get('misc', 'legend').copy()
    if settings.has_option('misc', 'plegend'):
        leg.update(settings.get('misc', 'plegend'))

    try:
        kwargs['best_peptides'] = settings.get('scoring', 'best peptides')
    except:
        kwargs['best_peptides'] = False

    best_res_raw, best_res = utils.multimap(n, func, peps, global_data=global_data, **kwargs)

    t2s_global = {}
    for global_data_local in global_data:
        t2s_global.update(global_data_local['t2s'])

    for spec_t, v in best_res_raw.items():
        peptide, m, snp_label, score, st, c, info = v
        spec_results[spec_t]['spectrum'] = t2s_global[spec_t]
        info['pep_nm'] = m
        info['charge'] = c
        spec_results[spec_t]['top_scores'] = -score
        spec_results[spec_t]['sequences'] = peptide
        spec_results[spec_t]['info'] = info
        spec_results[spec_t]['snp_label'] = snp_label

    maxlen = settings.getint('search', 'peptide maximum length')
    dtype = np.dtype([('score', np.float64),
        ('seq', np.str_, maxlen + 2), ('note', np.str_, 1),
        ('charge', np.int8), ('info', np.object_), ('sumI', np.float64), ('fragmentMT', np.float64), ('snp_label', np.str_, 15), ('nextscore_std', np.float64)])
    for spec_name, val in spec_results.items():
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
