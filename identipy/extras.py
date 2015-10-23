from scipy.stats import percentileofscore, scoreatpercentile
from pyteomics import achrom, auxiliary, parser, mass
from main import *
from scoring import get_fragment_mass_tol
import numpy as np
from utils import get_info, get_aa_mass, get_output, get_enzyme
try:
    from pyteomics import cmass
except ImportError:
    cmass = mass

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


def get_cutoff(results, settings, FDR=1):
    """A function for e-value threshold calculation"""
    target_evalues, decoy_evalues = [], []
    for res in get_output(results, settings):
        for e, (_, _, note, _, _, _) in zip(res['e-values'], res['candidates']):
            if note == 't':
                target_evalues.append(float(e))
            elif note == 'd':
                decoy_evalues.append(float(e))
    target_evalues = np.array(target_evalues)
    decoy_evalues = np.array(decoy_evalues)
    target_evalues.sort()
    best_cut_evalue = None
    for cut_evalue in target_evalues:
        counter_target = target_evalues[target_evalues <= cut_evalue].size
        counter_decoy = decoy_evalues[decoy_evalues <= cut_evalue].size
        if counter_target and (float(counter_decoy) / float(counter_target)) * 100 <= float(FDR):
            best_cut_evalue = cut_evalue
    if not best_cut_evalue:
        return 1e6
    return best_cut_evalue


def optimization(fname, settings):
    settings = copy(settings)
    settings.set('misc', 'first stage', '')
    settings.set('scoring', 'e-values for candidates', 1)

    results = []
    for res in process_file(fname, settings):
        results.append(res)
    print 'Results before optimization:', len(results)
    cutoff = get_cutoff(results, settings, FDR=1)
    print 'E-value cutoff:', cutoff

    functions = ['rt_filtering', 'precursor_mass_optimization', 'fragment_mass_optimization',
            'charge_optimization', 'missed_cleavages_optimization']
    for func in functions:
        settings = eval('%s(results, settings, cutoff)' % (func, ))
    return settings


def charge_optimization(results, settings, cutoff):
    settings = copy(settings)
    formula = "get_info(res['spectrum'], res, settings)[1]"
    chargestates = get_values(formula, results, settings, cutoff)
    mincharge = min(chargestates)
    maxcharge = max(chargestates)
    for ch in sorted(set(chargestates)):
        if float(chargestates[chargestates < ch].size) / chargestates.size < 0.05:
            mincharge = ch
    for ch in sorted(set(chargestates))[::-1]:
        if float(chargestates[chargestates > ch].size) / chargestates.size < 0.05:
            maxcharge = ch
    print 'NEW charges = %s:%s' % (mincharge, maxcharge)
    settings.set('search', 'maximum charge', maxcharge)
    settings.set('search', 'minimum charge', mincharge)
    return settings

def precursor_mass_optimization(results, settings, cutoff):
    settings = copy(settings)
    formula = """(cmass.fast_mass(str(seq), charge=0, aa_mass=get_aa_mass(settings)) - get_info(res['spectrum'], res, settings)[0]) / get_info(res['spectrum'], res, settings)[0] * 1e6"""
    massdif = get_values(formula, results, settings, cutoff)

    best_par_mt_l = min(massdif[massdif > scoreatpercentile(massdif, 0.5)])
    best_par_mt_r = max(massdif[massdif < scoreatpercentile(massdif, 99.5)])
    print 'NEW PARENT MASS TOLERANCE = %s:%s' % (best_par_mt_l, best_par_mt_r)
    settings.set('search', 'precursor accuracy left', -(best_par_mt_l))
    settings.set('search', 'precursor accuracy right', best_par_mt_r)
    return settings


def get_values(formula, results, settings, cutoff):
    values = np.array([])
    for res in get_output(results, settings):
        for e, (_, seq, note, _, _, _) in zip(res['e-values'], res['candidates']):
            e, seq, note, RT, _ = float(e), seq, note, utils.get_RT(res['spectrum']), res['spectrum']
            if note == 'd' or float(e) > cutoff:
                pass
            else:
                toadd = eval(formula)
                if isinstance(toadd, list):
                    if values.size:
                        values = np.vstack((values, toadd))
                    else:
                        values = np.array(toadd)
                else:
                    values = np.append(values, toadd)
    return values


def missed_cleavages_optimization(results, settings, cutoff):
    settings = copy(settings)
    formula = """len(parser.cleave(seq, get_enzyme(str(settings.get('search', 'enzyme'))), 0)) - 1"""

    missedcleavages = get_values(formula, results, settings, cutoff)
    best_missedcleavages = max(missedcleavages)
    for mc in sorted(set(missedcleavages))[::-1]:
        if float(missedcleavages[missedcleavages > mc].size) / missedcleavages.size < 0.05:
            best_missedcleavages = mc
    print 'NEW miscleavages = %s' % (best_missedcleavages, )
    settings.set('search', 'miscleavages', best_missedcleavages)
    return settings


def fragment_mass_optimization(results, settings, cutoff):
    settings = copy(settings)
    formula = """get_fragment_mass_tol(res['spectrum'], seq, settings)['fmt']"""
    fragmassdif = get_values(formula, results, settings, cutoff)
    step = FDbinSize(fragmassdif)
    lside, rside = 0, 1
    mt_h, _ = np.histogram(fragmassdif, bins=np.arange(lside, rside, step))

    for idx, mt in enumerate(mt_h):
        if mt == 0:
            mt_h = mt_h[:idx]
            break
#   threshold = mt_h.size * step
#   fragmassdif = fragmassdif[fragmassdif <= threshold]
    best_frag_mt = max(fragmassdif[fragmassdif < scoreatpercentile(fragmassdif, 97.5)])

    print 'NEW FRAGMENT MASS TOLERANCE = %s' % (best_frag_mt, )
    settings.set('search', 'product accuracy', best_frag_mt)
    return settings


def rt_filtering(results, settings, cutoff):
    settings = copy(settings)
    formula = "[float(RT), seq]"
    RTexp, seqs = zip(*get_values(formula, results, settings, cutoff))
    seqs = [list(s) for s in seqs] # FIXME: add terminal groups
    RTexp = [float(x) for x in RTexp]
    print len(RTexp), 'top PSMs with 1% FDR'
    RC_def = achrom.RCs_gilar_rp
    xdict = {}
    for key, val in RC_def['aa'].items():
        xdict[key] = [val, None]
    RC_dict = achrom.get_RCs_vary_lcp(seqs, RTexp)
    RC_dict_new = dict()
    for key, val in RC_dict['aa'].items():
        xdict.setdefault(key, [val, None])[1] = val
    a, b, _, _ = auxiliary.linear_regression([x[0] for x in xdict.values() if x[1] != None], [x[1] for x in xdict.values() if x[1] != None])
    for key, x in xdict.items():
        if x[1] == None:
            x[1] = x[0] * a + b
        RC_dict_new[key] = x[1]
    RC_dict['aa'] = RC_dict_new
    print RC_dict
    rtexp = np.array([np.mean(x) for x in RTexp])
    rttheor = np.array([achrom.calculate_RT(pep, RC_dict, raise_no_mod=False)
        for pep in seqs])
    deltaRT = rtexp - rttheor
    print auxiliary.linear_regression(rtexp, rttheor)

    def condition(spectrum, cand, _):
        return 1 <= percentileofscore(deltaRT, utils.get_RT(spectrum)
                        - achrom.calculate_RT(list(cand), RC_dict)
                    ) <= 99

    settings.set('scoring', 'condition', condition)
    return settings
