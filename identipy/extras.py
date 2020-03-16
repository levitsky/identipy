from scipy.stats import percentileofscore, scoreatpercentile
from scipy.optimize import curve_fit
from pyteomics import achrom, auxiliary as aux, parser, mass
from collections import Counter, defaultdict
from .main import *
from .scoring import get_fragment_mass_tol, get_fragment_mass_tol_ppm
import logging
logger = logging.getLogger(__name__)
import numpy as np
from .utils import get_info, get_aa_mass, get_enzyme, calculate_RT, get_title
try:
    from pyteomics import cmass
except ImportError:
    cmass = mass
from scipy.stats import rankdata
from copy import deepcopy
from scipy.optimize import curve_fit

def FDbinSize(X):
    """Calculates the Freedman-Diaconis bin size for
    a data set for use in making a histogram
    Arguments:
    X:  1D Data set
    Returns:
    h:  F-D bin size
    """
    X = np.sort(X)
    upperQuartile = scoreatpercentile(X, 75)
    lowerQuartile = scoreatpercentile(X, 25)
    IQR = upperQuartile - lowerQuartile
    h = 2. * IQR / len(X) ** (1. / 3.)
    return h

def get_peptides_subset(results):
    tmp_dict = dict()

    massdif = np.array([res['candidates'][0][4]['mzdiff']['ppm'] for res in results])

    for result in results:
        r_spectrum = get_title(result['spectrum'])
        r_sequence = str(result['candidates'][0][1])
        r_mass_diff_abs = abs(result['candidates'][0][4]['mzdiff']['ppm'])
        if r_sequence not in tmp_dict or r_mass_diff_abs < tmp_dict[r_sequence][1]:
            tmp_dict[r_sequence] = (r_spectrum, r_mass_diff_abs)
            # print(r_spectrum)

    new_results = []
    for result in results:
        r_spectrum = get_title(result['spectrum'])
        r_sequence = str(result['candidates'][0][1])
        if r_spectrum == tmp_dict[r_sequence][0]:
            new_results.append(result)
    return new_results

def get_subset(results, settings, fdr=0.01):
    """Filter results to given FDR using top 1 candidates"""
    subset = aux.filter(results, key=lambda x: x['e-values'][0],
            is_decoy = lambda x: x['candidates'][0][2] == 'd',
            fdr=fdr)
    return subset

def optimization(fname, settings):
    settings = settings.copy()
    settings.set('misc', 'first stage', '')
    efc = settings.get('scoring', 'e-values for candidates')
    settings.set('scoring', 'e-values for candidates', 1)
    left = settings.getfloat('search', 'precursor accuracy left')
    right = settings.getfloat('search', 'precursor accuracy right')
    wide = settings.getboolean('optimization', 'increase precursor mass tolerance')
    if settings.get('search', 'precursor accuracy unit') != 'ppm':
        left *= 1000
        right *= 1000
    if left < 100 and wide:
        settings.set('search', 'precursor accuracy left', 100)
    if right < 100 and wide:
        settings.set('search', 'precursor accuracy right', 100)
    # settings.set('search', 'precursor accuracy unit', 'ppm')
    results = list(process_file(fname, settings, initial_run=False))
    filtered = get_subset(results, settings, fdr=0.01)
    filtered = get_peptides_subset(filtered)
    logger.info('%s PSMs with 1%% FDR.', len(filtered))
    if len(filtered) < 250:
        if len(filtered) < 250:
            logger.warning('OPTIMIZATION ABORTED')
            return settings
        else:
            functions = [precursor_mass_optimization, fragment_mass_optimization,
                    missed_cleavages_optimization]
    else:
        functions = [
                rt_filtering,
                # precursor_mass_optimization,
                fragment_mass_optimization,
#               missed_cleavages_optimization
                ]
    for func in functions:
        # settings = func(filtered, settings, get_subset(results, settings, fdr=100.0))
        settings = func(filtered, settings, [x for x in results if x['candidates'][0][2] == 'd'])
    settings.set('scoring', 'e-values for candidates', efc)
    settings.set('scoring', 'best peptides', [str(res['candidates'][0][1]) for res in results])
    return settings


def charge_optimization(results, settings):
    settings = settings.copy()
    chargestates = np.array([get_info(res['spectrum'], res, settings)[1] for res in results])
    mincharge = chargestates.min()
    maxcharge = chargestates.max()
    
    for ch in range(mincharge, maxcharge+1):
        if float(chargestates[chargestates < ch].size) / chargestates.size < 0.01:
            mincharge = ch
    for ch in range(maxcharge, mincharge-1, -1):
        if float(chargestates[chargestates > ch].size) / chargestates.size < 0.01:
            maxcharge = ch
    logger.info('NEW charges = %s:%s', mincharge, maxcharge)
    settings.set('search', 'maximum charge', maxcharge)
    settings.set('search', 'minimum charge', mincharge)
    return settings

def calibrate_mass(bwidth, mass_left, mass_right, true_md):
    bbins = np.arange(-mass_left, mass_right, bwidth)
    H1, b1 = np.histogram(true_md, bins=bbins)
    b1 = b1 + bwidth
    b1 = b1[:-1]

    popt, pcov = curve_fit(noisygaus, b1, H1, p0=[1, np.median(true_md), 1, 1])
    mass_shift, mass_sigma = popt[1], np.abs(popt[2])
    return mass_shift, mass_sigma, pcov[0][0]

def noisygaus(x, a, x0, sigma, b):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + b

def precursor_mass_optimization(results, settings, unf):
    settings_nopime = settings.copy()
    settings_nopime.set('search', 'precursor isotope mass error', '0')
    settings_nopime.set('search', 'shifts', '0')
#   results = get_output(results, settings_nopime)

    settings = settings.copy()
    mass_left = settings.getfloat('search', 'precursor accuracy left')
    mass_right = settings.getfloat('search', 'precursor accuracy right')
    massdif = np.array([res['candidates'][0][4]['mzdiff']['ppm'] for res in results])
    massdif = massdif[(massdif > -mass_left) & (massdif < mass_right)]
    if settings.get('search', 'precursor accuracy unit') != 'ppm':
        mass_left = mass_left * 1e6 / 400
        mass_right = mass_right * 1e6 / 400
    logger.info('mass_left, mass_right: %s, %s', mass_left, mass_right)
    try:
        mass_shift, mass_sigma, covvalue = calibrate_mass(0.1, mass_left, mass_right, massdif)
        if np.isinf(covvalue):
            mass_shift, mass_sigma, covvalue = calibrate_mass(0.01, mass_left, mass_right, massdif)
        logger.info('%s, %s -> %s +- 8 * %s; %s', mass_left, mass_right, mass_shift, mass_sigma, covvalue)
        best_par_mt_l = mass_shift - 8 * mass_sigma
        best_par_mt_r = mass_shift + 8 * mass_sigma
        logger.info('SMART MASS TOLERANCE = %s:%s', best_par_mt_l, best_par_mt_r)
    except RuntimeError:
        error = True
    else:
        error = False
    if not error and np.isinf(covvalue):
        error = True
        logger.warning('Double error when fitting precursor errors: %s', massdif)
    print(percentileofscore(massdif, best_par_mt_r) - percentileofscore(massdif, best_par_mt_l), '!!!')
    if error or (percentileofscore(massdif, best_par_mt_r) - percentileofscore(massdif, best_par_mt_l) < 95):
        best_par_mt_l = scoreatpercentile(massdif, 0.1)
        best_par_mt_r = scoreatpercentile(massdif, 99.9)
        logger.warning('Percentage sanity check FAILED. Falling back on percentage boundaries')
    else:
        best_par_mt_l = max(best_par_mt_l, scoreatpercentile(massdif, 0.1))
        best_par_mt_r = min(best_par_mt_r, scoreatpercentile(massdif, 99.9))

    best_par_mt_l = -10
    best_par_mt_r = 10
    logger.info('NEW PARENT MASS TOLERANCE = %s:%s', best_par_mt_l, best_par_mt_r)
    settings.set('search', 'precursor accuracy left', -best_par_mt_l)
    settings.set('search', 'precursor accuracy right', best_par_mt_r)
    settings.set('search', 'precursor accuracy unit', 'ppm')
    return settings

def missed_cleavages_optimization(results, settings, unf):
    settings = settings.copy()
    missedcleavages = np.array([parser.num_sites(str(res['candidates'][0][1]), get_enzyme(str(settings.get('search', 'enzyme'))))
        for res in results])
    best_missedcleavages = missedcleavages.max()
    for mc in range(best_missedcleavages, -1, -1):
        if float(missedcleavages[missedcleavages > mc].size) / missedcleavages.size < 0.002:
            best_missedcleavages = mc
    logger.info('NEW miscleavages = %s', best_missedcleavages)
    settings.set('search', 'number of missed cleavages', best_missedcleavages)
    return settings

def fragment_mass_optimization(results, settings, results_unf):
    settings = settings.copy()
    fragmassdif = []
    I_all = []

    maxcharge = settings.getint('search', 'maximum charge')
    mincharge = settings.getint('search', 'minimum charge')

    rank_map = dict()
    ttl_cnt_d = dict()
    for ch in range(mincharge, maxcharge+1, 1):
        ttl_cnt_d[ch] = 0
        rank_map[ch] = dict()
        for i in range(1, 51, 1):
            rank_map[ch][i] = defaultdict(float)
    # ttl_cnt = 0

    # rank_map_unf = dict()
    # for i in range(1, 51, 1):
    #     rank_map_unf[i] = defaultdict(float)
    ttl_cnt_unf = 0

    # rank_map = dict()
    # for i in range(1, 51, 1):
    #     rank_map[i] = defaultdict(list)

    bions_map = dict()
    yions_map = dict()
    bions_map_all = dict()
    yions_map_all = dict()
    fragmassdif = []
    # fragmassdif_Da = []
    # logger.info(len(results))
    # logger.info(len(results_unf))
    if settings.has_option('misc', 'aa_mass'):
        aa_mass = settings.get('misc', 'aa_mass')
    else:
        aa_mass = get_aa_mass(settings)

    for res in results:

        neutral_mass, charge_state, RT = get_info(res['spectrum'], res, settings, aa_mass)
        ttl_cnt_d[charge_state] += 1
        p_len = len(str(res['candidates'][0][1]))
        tres = get_fragment_mass_tol(res['spectrum'], str(res['candidates'][0][1]), settings, charge_state)
        fragmassdif.extend(tres['fmt'])
        I_all.extend(tres['iall'])
        bion_cur = tres['bions']
        yion_cur = tres['yions']
        bion_curI = tres['bionsI']
        yion_curI = tres['yionsI']
        all_ion_cur = tres['allions']
        all_ion_curI = tres['allionsI']
        # print(bion_curI)

    #     # for ion in all_ion_cur:
    #     #     idx2 = 0
    #     #     for idx, b_cur in enumerate(all_ion_cur[ion]):
    #     #         ion_name = str(ion)+str(idx)
    #     #         if b_cur:
    #     #             i_rank = all_ion_curI[ion][idx2]
    #     #             rank_map[p_len][ion_name].append(i_rank)
    #     #             idx2 += 1

    #     # for ion in all_ion_cur:
    #     #     idx2 = 0
    #     #     for idx, b_cur in enumerate(all_ion_cur[ion]):
    #     #         ion_name = str(ion)+str(idx)
    #     #         if b_cur:
    #     #             i_rank = all_ion_curI[ion][idx2]
    #     #             if i_rank in rank_map[ion_name]:
    #     #                 rank_map[ion_name][i_rank] += 1.0
    #     #             else:
    #     #                 rank_map[ion_name][i_rank] = 1.0
    #     #             idx2 += 1
    #     #         else:
    #     #             if 'u' in rank_map[ion_name]:
    #     #                 rank_map[ion_name]['u'] += 1.0
    #     #             else:
    #     #                 rank_map[ion_name]['u'] = 1.0

    #     for ion in all_ion_cur:
    #         idx2 = 0
    #         for idx, b_cur in enumerate(all_ion_cur[ion]):
    #             if b_cur:
    #                 i_rank = all_ion_curI[ion][idx2]
    #                 rank_map[charge_state][i_rank][ion] += 1
    #                 idx2 += 1
    #     # for ion in all_ion_cur:
    #     #     idx2 = 0
    #     #     for idx, b_cur in enumerate(all_ion_cur[ion]):
    #     #         if b_cur:
    #     #             i_rank = all_ion_curI[ion][idx2]
    #     #             rank_map[i_rank][str(ion)+str(idx)] += 1
    #     #             idx2 += 1
    #     # for ion in all_ion_cur:
    #     #     idx2 = 0
    #     #     for idx, b_cur in enumerate(all_ion_cur[ion]):
    #     #         if b_cur:
    #     #             i_rank = all_ion_curI[ion][idx2]
    #     #             rank_map[i_rank][str(ion)+str(idx)] += 1
    #     #             idx2 += 1
    #     # idx2 = 0
    #     # for idx, y_cur in enumerate(yion_cur):
    #     #     if y_cur:
    #     #         i_rank = yion_curI[idx2]
    #     #         rank_map[i_rank]['y'+str(idx)] += 1

        

    #     # logger.info(bion_cur)
    #     if p_len not in bions_map:
    #         bions_map[p_len] = [0]*p_len
    #         yions_map[p_len] = [0]*p_len
    #     for idx, b_cur in enumerate(bion_cur):
    #         bions_map[p_len][idx] += b_cur
    #     for idx, y_cur in enumerate(yion_cur):
    #         yions_map[p_len][idx] += y_cur
    #     # logger.info(bion_cur)
    #     # logger.info(bion_curI)
    #     # logger.info(p_len)
    #     # logger.info(len(bion_cur))
    #     # if p_len not in bions_map:
    #     #     bions_map[p_len] = [[] for i in range(p_len-1)]
    #     #     yions_map[p_len] = [[] for i in range(p_len-1)]
    #     # idx2 = 0
    #     # for idx, b_cur in enumerate(bion_cur):
    #     #     if b_cur:
    #     #         bions_map[p_len][idx].append(bion_curI[idx2])
    #     #         idx2 += 1
    #     #     else:
    #     #         bions_map[p_len][idx].append(0)
    #     # idx2 = 0
    #     # for idx, y_cur in enumerate(yion_cur):
    #     #     if y_cur:
    #     #         yions_map[p_len][idx].append(yion_curI[idx2])
    #     #         idx2 += 1
    #     #     else:
    #     #         yions_map[p_len][idx].append(0)
    #     # for idx, (b_cur, b_curI) in enumerate(zip(bion_cur, bion_curI)):
    #     #     if b_cur:
    #     #         bions_map[p_len][idx].append(b_curI)
    #     # for idx, (y_cur, y_curI) in enumerate(zip(yion_cur, yion_curI)):
    #     #     if y_cur:
    #     #         yions_map[p_len][idx].append(y_curI)
    # # for k in rank_map:
    # #     for kk in rank_map[k]:
    # #         rank_map[k][kk] = rank_map[k][kk] / ttl_cnt
    # # for k in rank_map:
    # #     rank_map[k]['u'] = max(rank_map[k].values())

    # # for k in rank_map:
    # #     rank_map[k]['u'] = ttl_cnt - sum(rank_map[k].values())
    # #     for kk in rank_map[k]:
    # #         rank_map[k][kk] = rank_map[k][kk] / ttl_cnt

    # # for k in rank_map:
    # #     rank_map[k]['u'] = np.median(rank_map[k].values())
    # #     # rank_map[k]['u'] = min(rank_map[k].values())
    # #     koef = sum(rank_map[k].values())
    # #     for kk in rank_map[k]:
    # #         rank_map[k][kk] = rank_map[k][kk] / koef

    # for ch in rank_map:
    #     for k in rank_map[ch]:
    #         for kk in rank_map[ch][k]:
    #             rank_map[ch][k][kk] = rank_map[ch][k][kk] / ttl_cnt_d[ch]

    # try:
    #     import cPickle as pickle
    # except ImportError:
    #     import pickle
    # filenamep = '/home/mark/2020_poster_Denmark/rank_map.pickle'
    # with open(filenamep, 'wb') as output:
    #     pickle.dump(rank_map, output)

    # for ch in list(rank_map.keys()):
    #     all_vals = []
    #     for k in rank_map[ch]:
    #         all_vals.extend(list(rank_map[ch][k].values()))
    #     u_val = scoreatpercentile(all_vals, 1)
    #     if len(all_vals) == 0:
    #         del rank_map[ch]
    #     else:
    #         rank_map[ch]['u'] = u_val
    # for ch in rank_map:
    #     for k in rank_map[ch]:
    #         if k != 'u':
    #             for kk in list(rank_map[ch][k].keys()):
    #                 rank_map[ch][k][kk] = np.log(float(rank_map[ch][k][kk])/rank_map[ch]['u'])
    #     all_vals = []
    #     for key, val in rank_map[ch].items():
    #         if key != 'u':
    #             all_vals.extend(list(val.values()))
    # #     print(all_vals, ch)
    #     rank_map[ch]['m'] = min(all_vals)
    # for ch in range(mincharge, maxcharge+1, 1):
    #     if ch not in rank_map:
    #         try:
    #             rank_map[ch] = rank_map[ch-1]
    #         except:
    #             rank_map[ch] = rank_map[ch+1]



    # try:
    #     import cPickle as pickle
    # except ImportError:
    #     import pickle
    # # # print(rank_map[1])
    # # print(rank_map[2][1])
    # # print(rank_map[3][1])
    # # # print(rank_map[4])
    # # print(rank_map)

    # filenamep = '/home/mark/2020_poster_Denmark/rank_map.pickle'
    # with open(filenamep, 'wb') as output:
    #     pickle.dump(rank_map, output)

    # for kk in list(rank_map.keys()):
    #     tmp = rank_map[kk]
    #     if len(tmp):
    #         for ion, val in list(tmp.items()):
    #             tmp[ion] = len(val)
    #         # mval = max(tmp.values())
    #         # for ion, val in list(tmp.items()):
    #         #     if tmp[ion] < float(mval)/5:
    #         #         del tmp[ion]
    #         rank_map[kk] = tmp
    #     print(kk, 'OK')

    # for kk in list(rank_map.keys()):
    #     tmp = rank_map[kk]
    #     ttl = 0
    #     for ion, val in list(tmp.items()):
    #         ttl += len(val)
    #         tmp[ion] = np.mean(val)
    #     tmp_keys = list(tmp.keys())
    #     tmp_vals = list(tmp[k] for k in tmp_keys)
    #     tmp_ranked_vals = tmp_vals#rankdata(tmp_vals, method='ordinal')
    #     # tmp_ranked_vals = rankdata(tmp_vals, method='ordinal')
    #     for k, v in zip(tmp_keys, tmp_ranked_vals):
    #         tmp[k] = v
    #     rank_map[kk] = tmp
    #     print(kk, 'OK', ttl)

    # for k in rank_map:
    #     koef = sum(rank_map[k].values())
    #     for kk in rank_map[k]:
    #         rank_map[k][kk] = rank_map[k][kk] / koef

    # print(rank_map)


    # for res in results_unf:
    #     ttl_cnt_unf += 1
    #     tres = get_fragment_mass_tol(res['spectrum'], str(res['candidates'][0][1]), settings)
    #     all_ion_cur = tres['allions']
    #     all_ion_curI = tres['allionsI']

    #     # for ion in all_ion_cur:
    #     #     idx2 = 0
    #     #     for idx, b_cur in enumerate(all_ion_cur[ion]):
    #     #         if b_cur:
    #     #             i_rank = all_ion_curI[ion][idx2]
    #     #             rank_map_unf[i_rank][str(ion)+str(idx)] += 1
    #     #             idx2 += 1

    #     for ion in all_ion_cur:
    #         idx2 = 0
    #         for idx, b_cur in enumerate(all_ion_cur[ion]):
    #             ion_name = str(ion)+str(idx)
    #             if b_cur:
    #                 i_rank = all_ion_curI[ion][idx2]
    #                 if i_rank in rank_map_unf[ion_name]:
    #                     rank_map_unf[ion_name][i_rank] += 1.0
    #                 else:
    #                     rank_map_unf[ion_name][i_rank] = 1.0
    #                 idx2 += 1
    #             else:
    #                 if 'u' in rank_map_unf[ion_name]:
    #                     rank_map_unf[ion_name]['u'] += 1.0
    #                 else:
    #                     rank_map_unf[ion_name]['u'] = 1.0

    # filenamep = '/home/mark/2020_poster_Denmark/rank_map_unf.pickle'
    # with open(filenamep, 'wb') as output:
    #     pickle.dump(rank_map_unf, output)

    # for k in rank_map_unf:
    #     koef = sum(rank_map_unf[k].values())
    #     for kk in rank_map_unf[k]:
    #         rank_map_unf[k][kk] = rank_map_unf[k][kk] / koef

    # print(rank_map)
    # # for k in rank_map_unf:
    # #     for kk in rank_map_unf[k]:
    # #         rank_map_unf[k][kk] = rank_map_unf[k][kk] / ttl_cnt_unf
    # # for k in rank_map_unf:
    # #     rank_map_unf[k]['u'] = max(rank_map_unf[k].values())


    # for res in results_unf:
    #     p_len = len(str(res['candidates'][0][1]))
    #     bion_cur = get_fragment_mass_tol(res['spectrum'], str(res['candidates'][0][1]), settings)['bions']
    #     yion_cur = get_fragment_mass_tol(res['spectrum'], str(res['candidates'][0][1]), settings)['yions']
    #     # logger.info(bion_cur)
    #     if p_len not in bions_map_all:
    #         bions_map_all[p_len] = [0]*p_len
    #         yions_map_all[p_len] = [0]*p_len
    #     for idx, b_cur in enumerate(bion_cur):
    #         bions_map_all[p_len][idx] += b_cur
    #     for idx, y_cur in enumerate(yion_cur):
    #         yions_map_all[p_len][idx] += y_cur
    # # logger.info(yions_map)


    # for k in bions_map.keys():
    #     for idx, kk in enumerate(bions_map[k]):
    #         bions_map[k][idx] = np.mean(kk) 
    #     bions_map[k] = np.array(bions_map[k], dtype=float)
    #     bions_map[k] = np.nan_to_num(bions_map[k])

    # for k in yions_map.keys():
    #     for idx, kk in enumerate(yions_map[k]):
    #         yions_map[k][idx] = np.mean(kk) 
    #     yions_map[k] = np.array(yions_map[k], dtype=float)
    #     yions_map[k] = np.nan_to_num(yions_map[k])


    # for k in bions_map.keys():
    #     bions_map[k] = np.array(bions_map[k], dtype=float)
    #     koef = bions_map[k].sum()
    #     bions_map[k] = bions_map[k]/koef
    #     bions_map[k] = np.nan_to_num(bions_map[k])

    #     if k in bions_map_all:
    #         bions_map_all[k] = np.array(bions_map_all[k], dtype=float)
    #         koef = bions_map_all[k].sum()
    #         bions_map_all[k] = bions_map_all[k]/koef
    #         bions_map_all[k] = np.nan_to_num(bions_map_all[k])
    #     else:
    #         bions_map_all[k] = np.array([0]*k)
    #         yions_map_all[k] = np.array([0]*k)

    #     # bions_map[k] = bions_map[k] / bions_map_all[k]
    #     # bions_map[k] = np.nan_to_num(bions_map[k])

    #     bions_map[-k] = bions_map_all[k]
        
    # for k in yions_map.keys():
    #     yions_map[k] = np.array(yions_map[k], dtype=float)
    #     koef = yions_map[k].sum()
    #     yions_map[k] = yions_map[k]/koef
    #     yions_map[k] = np.nan_to_num(yions_map[k])

    #     yions_map_all[k] = np.array(yions_map_all[k], dtype=float)
    #     koef = yions_map_all[k].sum()
    #     yions_map_all[k] = yions_map_all[k]/koef
    #     yions_map_all[k] = np.nan_to_num(yions_map_all[k])

    #     # yions_map[k] = yions_map[k] / yions_map_all[k]
    #     # yions_map[k] = np.nan_to_num(yions_map[k])

    #     yions_map[-k] = yions_map_all[k]
    # logger.info(yions_map)
    fragmassdif = np.array(fragmassdif)
    # fragmassdif_Da = np.array(fragmassdif_Da)

    # try:
    #     import cPickle as pickle
    # except ImportError:
    #     import pickle
    # filenamep = '/home/mark/2020_poster_Denmark/frag.pickle'
    # with open(filenamep, 'wb') as output:
    #     pickle.dump(fragmassdif, output)

    # filenamep = '/home/mark/2020_poster_Denmark/iall.pickle'
    # with open(filenamep, 'wb') as output:
    #     pickle.dump(I_all, output)


    # print(len(I_all), len(fragmassdif))
    # return        

    best_frag_mt = scoreatpercentile(fragmassdif, 68) * 4    
    # best_frag_mt = scoreatpercentile(fragmassdif, 99)    

    logger.info('NEW FRAGMENT MASS TOLERANCE ppm = %s', best_frag_mt)
    settings.set('search', 'product accuracy ppm', best_frag_mt)
    settings.set('search', 'product accuracy unit', 'ppm')

    # orig_acc = settings.getfloat('search', 'product accuracy')
    # settings.set('search', 'product accuracy', 0.5)

    try:
    # if 1:

        rank_map = {
            0: {},
            1: {},
        }
        ttl_cnt_d = {
            0: {},
            1: {},
        }
        
        for ch in range(mincharge, maxcharge+1, 1):
            ttl_cnt_d[0][ch] = 0
            ttl_cnt_d[1][ch] = 0
            rank_map[0][ch] = dict()
            rank_map[1][ch] = dict()
            for i in range(1, 51, 1):
                rank_map[0][ch][i] = defaultdict(float)
                rank_map[1][ch][i] = defaultdict(float)
        for res in results:
            # print(res['spectrum'])
            if 'params' in res['spectrum'] and 'isowidthdiff' in res['spectrum']['params'] and abs(float(res['spectrum']['params']['isowidthdiff'])) >= 0.1:
                chimeric = 1
            else:
                chimeric = 0
            neutral_mass, charge_state, RT = get_info(res['spectrum'], res, settings, aa_mass)
            ttl_cnt_d[chimeric][charge_state] += 1
            p_len = len(str(res['candidates'][0][1]))
            tres = get_fragment_mass_tol(res['spectrum'], str(res['candidates'][0][1]), settings, charge_state)
            # tres = get_fragment_mass_tol_ppm(res['spectrum'], str(res['candidates'][0][1]), settings, charge_state, acc_ppm=best_frag_mt)
            all_ion_cur = tres['allions']
            all_ion_curI = tres['allionsI']
            for ion in all_ion_cur:
                idx2 = 0
                for idx, b_cur in enumerate(all_ion_cur[ion]):
                    if b_cur:
                        i_rank = all_ion_curI[ion][idx2]
                        rank_map[chimeric][charge_state][i_rank][ion] += 1
                        idx2 += 1
        for chim in rank_map:
            for ch in rank_map[chim]:
                for k in rank_map[chim][ch]:
                    for kk in rank_map[chim][ch][k]:
                        rank_map[chim][ch][k][kk] = rank_map[chim][ch][k][kk] / ttl_cnt_d[chim][ch]

        try:
            import cPickle as pickle
        except ImportError:
            import pickle
        filenamep = '/home/mark/2020_poster_Denmark/IPGF1.pickle'
        with open(filenamep, 'wb') as output:
            pickle.dump(rank_map, output)

        # for ch in list(rank_map.keys()):
        #     all_vals = []
        #     for k in rank_map[ch]:
        #         rank_map[ch][k]['um'] = 1 - sum(rank_map[ch][k].values())
        #     #     all_vals.extend(list(rank_map[ch][k].values()))
        #     # u_val = scoreatpercentile(all_vals, 1)
        #     # if len(all_vals) == 0:
        #     #     del rank_map[ch]
        #     # else:
        #     #     rank_map[ch]['u'] = u_val
        #         # rank_map[ch]['u'] = 0.01
        # for ch in rank_map:
        #     for k in rank_map[ch]:
        #         if k != 'u':
        #             for kk in list(rank_map[ch][k].keys()):
        #                 if kk != 'um':
        #                     rank_map[ch][k][kk] = np.log(float(rank_map[ch][k][kk])/rank_map[ch][k]['um'])
        #                 # rank_map[ch][k][kk] = float(rank_map[ch][k][kk])
        #                 # rank_map[ch][k][kk] = max(0, np.log(float(rank_map[ch][k][kk])/rank_map[ch]['u']))
        #                 # rank_map[ch][k][kk] = np.log(float(rank_map[ch][k][kk])/rank_map[ch]['u'])
        #     # all_vals = []
        #     # for key, val in rank_map[ch].items():
        #     #     if key != 'u':
        #     #         all_vals.extend(list(val.values()))
        #     # rank_map[ch]['m'] = min(all_vals) - np.log(1./2.)

        for chim in list(rank_map.keys()):
            for ch in list(rank_map[chim].keys()):
                all_vals = []
                for k in rank_map[chim][ch]:
                    all_vals.extend(list(rank_map[chim][ch][k].values()))

            # noise_mean, noise_sigma, covvalue = calibrate_mass(0.01, 0, 1.0, all_vals)


                u_val = scoreatpercentile(all_vals, 1)
                if len(all_vals) == 0:
                    del rank_map[chim][ch]
                else:
                    # rank_map[ch]['u'] = np.median(all_vals) * 2
                    # rank_map[ch]['u'] = noise_mean + 3 * noise_sigma
                    rank_map[chim][ch]['u'] = u_val
                    # rank_map[ch]['u'] = 0.01
                    # rank_map[chim][ch]['u'] = np.std(all_vals)
                    # rank_map[chim][ch]['med'] = np.median(all_vals)
        for chim in rank_map:
            for ch in rank_map[chim]:
                for k in rank_map[chim][ch]:
                    if k != 'u' and k != 'med':
                        for kk in list(rank_map[chim][ch][k].keys()):
                            # rank_map[chim][ch][k][kk] = (float(rank_map[chim][ch][k][kk]) - rank_map[chim][ch]['med']) / rank_map[chim][ch]['u'] 
                            # rank_map[ch][k][kk] = max(0, np.log(float(rank_map[ch][k][kk])/rank_map[ch]['u'])) 
                            # min_val = min(rank_map[ch][k].values())
                            # if not min_val:
                            #     rank_map[ch][k][kk] = 0
                            # else:
                                # rank_map[ch][k][kk] = max(0, np.log(float(rank_map[ch][k][kk])/min_val))
                            rank_map[chim][ch][k][kk] = np.log(float(rank_map[chim][ch][k][kk])/rank_map[chim][ch]['u'])
                            # rank_map[chim][ch][k][kk] = max(0, np.log(float(rank_map[chim][ch][k][kk])/rank_map[chim][ch]['u']))
                all_vals = []
                for key, val in rank_map[chim][ch].items():
                    if key != 'u' and key != 'med':
                        all_vals.extend(list(val.values()))
                rank_map[chim][ch]['m'] = min(all_vals)# - np.log(1./4.)
        for ch in range(mincharge, maxcharge+1, 1):
            for chim in [0, 1]:
                if ch not in rank_map[chim]:
                    try:
                        rank_map[chim][ch] = rank_map[chim][ch-1]
                    except:
                        try:
                            rank_map[chim][ch] = rank_map[chim][ch+1]
                        except:
                            print('missing autofill for %d chim, %d charge' % (chim, ch))
        # rank_map[1] = rank_map[2]
        # rank_map[3] = rank_map[2]
        # rank_map[4] = rank_map[2]
        # rank_map[5] = rank_map[2]
        # print(rank_map[0][2][1])
        # try:
        #     print(rank_map[1][2][1])
        # except:
        #     pass
        # print(rank_map[4])
        # print(rank_map)


        # settings.set('search', 'product accuracy', orig_acc)


        try:
            import cPickle as pickle
        except ImportError:
            import pickle
        filenamep = '/home/mark/2020_poster_Denmark/IPGF1.pickle'
        with open(filenamep, 'wb') as output:
            pickle.dump(rank_map, output)


        settings.set('search', 'bions_map', bions_map)
        settings.set('search', 'yions_map', yions_map)
        settings.set('search', 'rank_map', deepcopy(rank_map))
        settings.set('search', 'rank_map_unf', rank_map)
        # logger.info(settings.get('search', 'bions_map'))


        rank_map = {
            0: {},
            1: {},
        }

        rank_mapw = {
            0: {},
            1: {},
        }

        ttl_cnt_d = {
            0: {},
            1: {},
        }


        ttl_cnt_ions = defaultdict(float)
        ttl_cnt_ions_all = defaultdict(float)
        ttl_res = 0
        
        for ch in range(mincharge, maxcharge+1, 1):
            ttl_cnt_d[0][ch] = 0
            ttl_cnt_d[1][ch] = 0
            rank_map[0][ch] = dict()
            rank_map[1][ch] = dict()
            rank_mapw[0][ch] = dict()
            rank_mapw[1][ch] = dict()
            for i in range(1, 51, 1):
                rank_map[0][ch][i] = defaultdict(float)
                rank_map[1][ch][i] = defaultdict(float)
                rank_mapw[0][ch][i] = defaultdict(float)
                rank_mapw[1][ch][i] = defaultdict(float)
        for res in results:
            ttl_res += 1
            # print(res['spectrum'])
            if 'params' in res['spectrum'] and 'isowidthdiff' in res['spectrum']['params'] and abs(float(res['spectrum']['params']['isowidthdiff'])) >= 0.1:
                chimeric = 1
            else:
                chimeric = 0
            neutral_mass, charge_state, RT = get_info(res['spectrum'], res, settings, aa_mass)
            ttl_cnt_d[chimeric][charge_state] += 1
            p_len = len(str(res['candidates'][0][1]))
            tres = get_fragment_mass_tol(res['spectrum'], str(res['candidates'][0][1]), settings, charge_state)
            # tres = get_fragment_mass_tol_ppm(res['spectrum'], str(res['candidates'][0][1]), settings, charge_state, acc_ppm=best_frag_mt)
            orig_acc = settings.getfloat('search', 'product accuracy')
            settings.set('search', 'product accuracy', 2 * orig_acc)
            tresw = get_fragment_mass_tol(res['spectrum'], str(res['candidates'][0][1]), settings, charge_state)
            settings.set('search', 'product accuracy', orig_acc)
            all_ion_cur = tres['allions']
            all_ion_curI = tres['allionsI']
            for ion in all_ion_cur:
                idx2 = 0
                for idx, b_cur in enumerate(all_ion_cur[ion]):
                    if b_cur:
                        i_rank = all_ion_curI[ion][idx2]
                        rank_map[chimeric][charge_state][i_rank][ion] += 1
                        idx2 += 1
                        ttl_cnt_ions[ion] += 1
                    ttl_cnt_ions_all[ion] += 1

            all_ion_curw = tresw['allions']
            all_ion_curIw = tresw['allionsI']
            for ion in all_ion_curw:
                idx2 = 0
                for idx, b_cur in enumerate(all_ion_curw[ion]):
                    if b_cur:
                        i_rank = all_ion_curIw[ion][idx2]
                        rank_mapw[chimeric][charge_state][i_rank][ion] += 1
                        idx2 += 1

        for k in list(ttl_cnt_ions.keys()):
            # ttl_cnt_ions[k] = ttl_cnt_ions[k] / ttl_res * 100
            ttl_cnt_ions[k] = ttl_cnt_ions[k]# / sum(vv for kk, vv in ttl_cnt_d[0].items() if max(1, (kk - 1)) >= k[-1]) * 100
            ttl_cnt_ions[k] = ttl_cnt_ions[k] / ttl_cnt_ions_all[k] * 100
        print(ttl_cnt_ions)


        try:
            import cPickle as pickle
        except ImportError:
            import pickle
        filenamep = '/home/mark/2020_poster_Denmark/ttl_cnt_ions_highres.pickle'
        with open(filenamep, 'wb') as output:
            pickle.dump(ttl_cnt_ions, output)

        massdif = list(ttl_cnt_ions.values())
        try:
            mass_shift, mass_sigma, covvalue = calibrate_mass(0.1, 0, 100, massdif)
            print(mass_shift, mass_sigma, covvalue)
            perc_threshold = mass_shift + 3 * mass_sigma
            print(perc_threshold)
        except:
            perc_threshold = 5
            print(perc_threshold)
        if perc_threshold < 0:
            perc_threshold = 5
            print(perc_threshold)

        allowed_ions = set()
        for k in ttl_cnt_ions:
            if ttl_cnt_ions[k] >= perc_threshold:
                allowed_ions.add(k)
        settings.set('search', 'allowed_ions', allowed_ions)


        for chim in list(rank_mapw.keys()):
            for ch in list(rank_mapw[chim].keys()):
                for k in rank_mapw[chim][ch]:
                    for kk in list(rank_mapw[chim][ch][k].keys()):
                        if kk not in allowed_ions or rank_mapw[chim][ch][k][kk] <= 5:
                            # print('low vals for %d %d %s %d' % (ch, k, kk, rank_mapw[chim][ch][k][kk]))
                            # rank_mapw[chim][ch][k][kk] = 0
                            # rank_map[chim][ch][k][kk] = 0
                            rank_mapw[chim][ch][k][kk] = 0
                            rank_map[chim][ch][k][kk] = 0
                        # else:
                            # print('high vals for %d %d %s %d' % (ch, k, kk, rank_mapw[chim][ch][k][kk]))

        for chim in rank_map:
            for ch in rank_map[chim]:
                for k in rank_map[chim][ch]:
                    for kk in rank_map[chim][ch][k]:
                        rank_map[chim][ch][k][kk] = rank_map[chim][ch][k][kk] / ttl_cnt_d[chim][ch]

        for chim in rank_mapw:
            for ch in rank_mapw[chim]:
                for k in rank_mapw[chim][ch]:
                    for kk in rank_mapw[chim][ch][k]:
                        rank_mapw[chim][ch][k][kk] = rank_mapw[chim][ch][k][kk] / ttl_cnt_d[chim][ch]

        # try:
        #     import cPickle as pickle
        # except ImportError:
        #     import pickle
        # filenamep = '/home/mark/2020_poster_Denmark/IPGF1.pickle'
        # with open(filenamep, 'wb') as output:
        #     pickle.dump(rank_map, output)

        # for ch in list(rank_map.keys()):
        #     all_vals = []
        #     for k in rank_map[ch]:
        #         rank_map[ch][k]['um'] = 1 - sum(rank_map[ch][k].values())
        #     #     all_vals.extend(list(rank_map[ch][k].values()))
        #     # u_val = scoreatpercentile(all_vals, 1)
        #     # if len(all_vals) == 0:
        #     #     del rank_map[ch]
        #     # else:
        #     #     rank_map[ch]['u'] = u_val
        #         # rank_map[ch]['u'] = 0.01
        # for ch in rank_map:
        #     for k in rank_map[ch]:
        #         if k != 'u':
        #             for kk in list(rank_map[ch][k].keys()):
        #                 if kk != 'um':
        #                     rank_map[ch][k][kk] = np.log(float(rank_map[ch][k][kk])/rank_map[ch][k]['um'])
        #                 # rank_map[ch][k][kk] = float(rank_map[ch][k][kk])
        #                 # rank_map[ch][k][kk] = max(0, np.log(float(rank_map[ch][k][kk])/rank_map[ch]['u']))
        #                 # rank_map[ch][k][kk] = np.log(float(rank_map[ch][k][kk])/rank_map[ch]['u'])
        #     # all_vals = []
        #     # for key, val in rank_map[ch].items():
        #     #     if key != 'u':
        #     #         all_vals.extend(list(val.values()))
        #     # rank_map[ch]['m'] = min(all_vals) - np.log(1./2.)


        for chim in list(rank_mapw.keys()):
            for ch in list(rank_mapw[chim].keys()):
                all_vals = []
                for k in rank_mapw[chim][ch]:
                    all_vals.extend(list(rank_mapw[chim][ch][k].values()))
                if len(all_vals) == 0:
                    del rank_mapw[chim][ch]
                    del rank_map[chim][ch]

        for chim in rank_mapw:
            for ch in rank_mapw[chim]:
                for k in rank_mapw[chim][ch]:
                    if k != 'u' and k != 'med':
                        for kk in list(rank_mapw[chim][ch][k].keys()):
                            # rank_map[chim][ch][k][kk] = (float(rank_map[chim][ch][k][kk]) - rank_map[chim][ch]['med']) / rank_map[chim][ch]['u'] 
                            # rank_map[ch][k][kk] = max(0, np.log(float(rank_map[ch][k][kk])/rank_map[ch]['u'])) 
                            # min_val = min(rank_map[ch][k].values())
                            # if not min_val:
                            #     rank_map[ch][k][kk] = 0
                            # else:
                                # rank_map[ch][k][kk] = max(0, np.log(float(rank_map[ch][k][kk])/min_val))
                            # rank_map[chim][ch][k][kk] = 1 - rank_map[chim][ch][k][kk]/rank_mapw[chim][ch][k][kk]

                            orig_v = rank_map[chim][ch][k].get(kk, 0)
                            ext_v = rank_mapw[chim][ch][k][kk]
                            false_v = float(ext_v - orig_v)# / 9
                            true_v = orig_v - false_v
                            # false_v = ext_v#float(ext_v - orig_v)# / 9
                            # true_v = orig_v# - false_v


                            if false_v and true_v:
                                rank_map[chim][ch][k][kk] = np.log(true_v / false_v)
                                # rank_map[chim][ch][k][kk] = -np.log(true_v / false_v)
                                # rank_map[chim][ch][k][kk] = np.sqrt(true_v / false_v)
                            else:
                                # rank_map[chim][ch][k][kk] = np.log(1e4)
                                # rank_map[chim][ch][k][kk] = np.sqrt(1e4)
                                if kk in rank_map[chim][ch][k]:
                                    del rank_map[chim][ch][k][kk]

                            # rank_map[chim][ch][k][kk] = -np.log(0.99999 - rank_map[chim][ch][k].get(kk, 0)/rank_mapw[chim][ch][k][kk])
                            # rank_map[chim][ch][k][kk] = rank_map[chim][ch][k].get(kk, 0)/rank_mapw[chim][ch][k][kk]


        for chim in list(rank_map.keys()):
            for ch in list(rank_map[chim].keys()):
                for k in list(rank_map[chim][ch].keys()):
                    all_vals = []
                    all_vals.extend(list(rank_map[chim][ch][k].values()))
                    if len(all_vals) == 0:
                        del rank_mapw[chim][ch][k]
                        del rank_map[chim][ch][k]

        for chim in list(rank_map.keys()):
            for ch in list(rank_map[chim].keys()):
                all_vals = []
                for k in rank_map[chim][ch]:
                    all_vals.extend(list(rank_map[chim][ch][k].values()))
                if len(all_vals) <= 5:
                    print('NOT OK %d %d' % (len(all_vals), ch))
                    del rank_mapw[chim][ch]
                    del rank_map[chim][ch]
                else:
                    print('OK %d %d' % (len(all_vals), ch))


        for chim in rank_map:
            for ch in rank_map[chim]:
                all_vals = []
                for key, val in rank_map[chim][ch].items():
                    if key != 'u' and key != 'med':
                        all_vals.extend(list(val.values()))
                rank_map[chim][ch]['m'] = min(all_vals)# - np.log(1./4.)#0#min(all_vals)# - np.log(1./4.)
                            
        for chim in list(rank_map.keys()):
            for ch in list(rank_map[chim].keys()):
                for k in list(range(1, max(z for z in rank_map[chim][ch] if type(z) != str)+1, 1))[::-1]:
                    if k not in rank_map[chim][ch]:
                        if k-1 in rank_map[chim][ch]:
                            rank_map[chim][ch][k] = rank_map[chim][ch][k-1]
                    if k in list(rank_map[chim][ch].keys()):
                        for kk in list(rank_map[chim][ch][1].keys()):
                            if kk not in rank_map[chim][ch][k]:
                                if k-1 in rank_map[chim][ch]:
                                    if kk in rank_map[chim][ch][k-1]:
                                        rank_map[chim][ch][k][kk] = rank_map[chim][ch][k-1][kk]

        for ch in range(mincharge, maxcharge+1, 1):
            for chim in [0, 1]:
                if ch not in rank_map[chim]:
                    try:
                        rank_map[chim][ch] = rank_map[chim][ch-1]
                    except:
                        try:
                            rank_map[chim][ch] = rank_map[chim][ch+1]
                        except:
                            print('missing autofill for %d chim, %d charge' % (chim, ch))
        # print(rank_map[0][2][1])
        # print(rank_map[0][2][5])

        try:
            import cPickle as pickle
        except ImportError:
            import pickle
        filenamep = '/home/mark/2020_poster_Denmark/IPGF2.pickle'
        with open(filenamep, 'wb') as output:
            pickle.dump(rank_map, output)

        settings.set('search', 'rank_map_unf', deepcopy(rank_map))

    except:
        print('IPGF scores deactivated!')
        settings.set('search', 'rank_map', False)
        settings.set('search', 'rank_map_unf', False)
        settings.set('search', 'bions_map', False)
        settings.set('search', 'yions_map', False)
        settings.set('search', 'allowed_ions', set())
        settings.set('search', 'use_allowed_ions', 0)


    return settings


def rt_filtering(results, settings, unf):
    settings = settings.copy()
    if settings.has_option('misc', 'legend'):
        legend = settings.get('misc', 'legend')
    else:
        legend = None
    RTexp, seqs = zip(*[(utils.get_RT(res['spectrum']), res['candidates'][0][1]) for res in results])
    if legend is not None:
        stdl = set(parser.std_labels)
        newseqs = []
        for s in seqs:
            if parser.fast_valid(s):
                newseqs.append(list(s))
            else:
                seq = []
                c, n = False, False
                for c in s:
                    if c in stdl:
                        seq.append(c)
                    else:
                        mod, res, term = legend[c]
                        if res == '-':
                            if term == '[':
                                seq.append(mod+'-')
                                n = True
                            else:
                                seq.append('-'+mod)
                                c = True
                        else:
                            seq.append(mod+res)
                    if not n: seq.append(parser.std_nterm)
                    if not c: seq.append(parser.std_cterm)
                newseqs.append(seq)
        seqs = newseqs
    RTexp = [float(x) for x in RTexp]
    if np.allclose(RTexp, 0):
        logger.warning('RT is missing. Skipping RT optimization.')
        return settings
    RC_def = achrom.RCs_gilar_rp
    xdict = {}
    for key, val in RC_def['aa'].items():
        xdict[key] = [val, None]
    RC_dict = utils.get_RCs_vary_lcp(seqs, RTexp)
    RC_dict_new = dict()
    for key, val in RC_dict['aa'].items():
        xdict.setdefault(key, [val, None])[1] = val
    a, b, _, _ = aux.linear_regression([x[0] for x in xdict.values() if x[1] != None], [x[1] for x in xdict.values() if x[1] != None])
    for key, x in xdict.items():
        if x[1] == None:
            x[1] = x[0] * a + b
        RC_dict_new[key] = x[1]
    if legend is not None:
        for k, v in legend.items():
            if len(k) == 1: continue
            if k[-1] in '[]':
                if k[-2] == '-':
                    kk = ('-' + k[1:-1]) if k[-1] == ']' else (k[:-1])
                else:
                    kk = k[:-1]
            elif len(k) > 1:
                kk = k
            logger.debug('%s -> %s', k, kk)
            if kk in RC_dict_new:
                RC_dict_new[v] = RC_dict_new[kk]
            else:
                if kk[-1].isupper():
                    kkk = kk[-1]
                elif kk[-1] == '-':
                    kkk = parser.std_nterm
                elif kk[0] == '-':
                    kkk = parser.std_cterm
                RC_dict_new[v] = RC_dict_new.get(kkk, 0)
                logger.info('No RC for %s, using %s or 0: %s', kk, kkk, RC_dict_new[v])


    RC_dict['aa'] = RC_dict_new

    logger.debug('RC dict: %s', RC_dict)
    rtexp = np.array([np.mean(x) for x in RTexp])
    rttheor = np.array([calculate_RT(pep, RC_dict, raise_no_mod=False)
        for pep in seqs])
    deltaRT = rtexp - rttheor
    logger.debug('Linear regression: %s', aux.linear_regression(rtexp, rttheor))
    best_RT_l = scoreatpercentile(deltaRT, 0.05)
    best_RT_r = scoreatpercentile(deltaRT, 99.95)

    def condition(spectrum, cand, _, stored_value=False):
        if not stored_value:
            stored_value = calculate_RT(cand, RC_dict)
        rtd = spectrum['RT'] - stored_value
        return best_RT_l <= rtd <= best_RT_r, stored_value
    settings.set('scoring', 'condition', condition)
    return settings


def calibrate_mass(bwidth, mass_left, mass_right, true_md):
    bbins = np.arange(-mass_left, mass_right, bwidth)
    H1, b1 = np.histogram(true_md, bins=bbins)
    b1 = b1 + bwidth
    b1 = b1[:-1]

    popt, pcov = curve_fit(noisygaus, b1, H1, p0=[1, np.median(true_md), 1, 1])
    mass_shift, mass_sigma = popt[1], np.abs(popt[2])
    return mass_shift, mass_sigma, pcov[0][0]

def noisygaus(x, a, x0, sigma, b):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + b
