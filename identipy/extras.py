from scipy.stats import percentileofscore, scoreatpercentile
from pyteomics import achrom, auxiliary
from main import *
from scoring import get_fragment_mass_tol
from numpy import mean, sort
from utils import get_info, get_aa_mass


def FDbinSize(X):
    """Calculates the Freedman-Diaconis bin size for
    a data set for use in making a histogram
    Arguments:
    X:  1D Data set
    Returns:
    h:  F-D bin size
    """
    X = sort(X)
    upperQuartile = scoreatpercentile(X, 75)
    lowerQuartile = scoreatpercentile(X, 25)
    IQR = upperQuartile - lowerQuartile
    h = 2. * IQR / len(X) ** (1. / 3.)
    return h


def get_cutoff(results, FDR=1):
    """A function for e-value threshold calculation"""
    target_evalues, decoy_evalues = np.array([]), np.array([])
    for res in results:
        for e, (_, _, note) in zip(res['e-values'], res['candidates']):
            if note == 't':
                target_evalues = np.append(target_evalues, float(e))
            elif note == 'd':
                decoy_evalues = np.append(decoy_evalues, float(e))
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


def smart_filtering(fname, settings):
    conditions = settings.get('smart filtering', 'conditions')
    cutoff = settings.getfloat('smart filtering', 'initial e-value cut-off')
    a, b = 1, 99
    # somehow derive the proper e-value cut-off
    # TODO
    settings = copy(settings)
    settings.set('misc', 'first stage', '')
    results = []
    for res in process_file(fname, settings):
        for e, (_, seq, note) in zip(res['e-values'], res['candidates']):
            results.append((float(e), seq, note, utils.get_RT(res['spectrum']), res['spectrum']))
    cutoff = get_cutoff(results, FDR=1)

    points = tuple([] for _ in conditions)
    for res in process_file(fname, settings):
        for e, (_, seq) in zip(res['e-values'], res['candidates']):
            if e > cutoff:
                break
            for p, (func, _, _) in zip(points, conditions):
                p.append(func(res['spectrum'], seq))

    def conjunction(spectrum, cand, _):
        return all(a <= percentileofscore(data, func(spectrum, cand)) <= b
                for data, (func, a, b) in zip(points, conditions))
    settings.set('scoring', 'condition', conjunction)
    return settings


def optimization(fname, settings):
    settings = copy(settings)
    settings.set('misc', 'first stage', '')
    settings.set('scoring', 'e-values for candidates', 1)

    results = []
    for res in process_file(fname, settings):
        results.append(res)
    cutoff = get_cutoff(results, FDR=1)
    print cutoff
    
    functions = ['rt_filtering', 'precursor_mass_optimization', 'fragment_mass_optimization']
    for func in functions:
        settings = eval('%s(results, settings, cutoff)' % (func, ))
    return settings

def precursor_mass_optimization(results, settings, cutoff):
    settings = copy(settings)
    massdif, seqs = np.array([]), []
    aa_mass = get_aa_mass(settings)
    j = 0
    for res in results:
        for e, (_, seq, note) in zip(res['e-values'], res['candidates']):
            e, seq, note, _, _ = float(e), seq, note, utils.get_RT(res['spectrum']), res['spectrum']
            j += 1
            if note == 'd' or float(e) > cutoff:
                pass
            else:
                neutral_mass, _, _ = get_info(res['spectrum'], res, settings)
                massdif = np.append(massdif, (neutral_mass - mass.fast_mass(seq, charge=0, aa_mass=aa_mass)) / neutral_mass * 1e6)
                seqs.append(str(seq))
    best_par_mt_l = min(massdif[massdif > scoreatpercentile(massdif, 0.5)])
    best_par_mt_r = max(massdif[massdif < scoreatpercentile(massdif, 99.5)])

    print 'NEW PARENT MASS TOLERANCE = %s:%s' % (best_par_mt_l, best_par_mt_r)
    settings.set('search', 'precursor accuracy value left', best_par_mt_l)
    settings.set('search', 'precursor accuracy value right', best_par_mt_r)
    return settings


def fragment_mass_optimization(results, settings, cutoff):
    settings = copy(settings)
    fragmassdif, seqs = np.array([]), []
    j = 0
    for res in results:
        for e, (_, seq, note) in zip(res['e-values'], res['candidates']):
            e, seq, note, _, _ = float(e), seq, note, utils.get_RT(res['spectrum']), res['spectrum']
            j += 1
            if note == 'd' or float(e) > cutoff:
                pass
            else:
                new_params = get_fragment_mass_tol(res['spectrum'], seq, settings)
                fragmassdif = np.append(fragmassdif, new_params['fmt'])
                seqs.append(str(seq))
    step = FDbinSize(fragmassdif)
    lside, rside = 0, 1
    mt_h, _ = np.histogram(fragmassdif, bins=np.arange(lside, rside, step))

    for idx, mt in enumerate(mt_h):
        if mt == 0:
            mt_h = mt_h[:idx]
            break
    threshold = mt_h.size * step
    fragmassdif = fragmassdif[fragmassdif <= threshold]
    best_frag_mt = max(fragmassdif[fragmassdif < scoreatpercentile(fragmassdif, 95)])

    print 'NEW FRAGMENT MASS TOLERANCE = %s' % (best_frag_mt, )
    settings.set('search', 'product accuracy', best_frag_mt)
    return settings
    
    
def rt_filtering(results, settings, cutoff):
    RTexp, seqs = [], []
    j = 0
    for res in results:
        for e, (_, seq, note) in zip(res['e-values'], res['candidates']):
            e, seq, note, RT, _ = float(e), seq, note, utils.get_RT(res['spectrum']), res['spectrum']
            j += 1
            if note == 'd' or float(e) > cutoff:
                pass
            else:
                RTexp.append(RT)
                seqs.append(str(seq))
    print len(RTexp), 'top PSMs with 1% FDR'
    RC_def = achrom.RCs_yoshida
    xdict = {}
    for key, val in RC_def['aa'].items():
        xdict[key] = [val, None]
    RC_dict = achrom.get_RCs_vary_lcp(seqs, RTexp)
    RC_dict_new = dict()
    for key, val in RC_dict['aa'].items():
        xdict[key][1] = val
    a, b, _, _ = auxiliary.linear_regression([x[0] for x in xdict.values() if x[1] != None], [x[1] for x in xdict.values() if x[1] != None])
    for key, x in xdict.items():
        if x[1] == None:
            x[1] = x[0] * a + b
        RC_dict_new[key] = x[1]
    RC_dict['aa'] = RC_dict_new
    deltaRT = [rtexp - achrom.calculate_RT(pep, RC_dict, raise_no_mod=False)
            for rtexp, pep in zip([mean(x) for x in RTexp], seqs)]

    def condition(spectrum, cand, _):
        return 1 <= percentileofscore(deltaRT, utils.get_RT(spectrum)
                        - achrom.calculate_RT(cand, RC_dict)
                    ) <= 99

    settings.set('scoring', 'condition', condition)
    return settings

def achrom_rt_filtering(fname, settings):
    settings = copy(settings)
    settings.set('misc', 'first stage', '')
    settings.set('scoring', 'e-values for candidates', 1)
    RTexp, massdif, fragstd, seqs = [], np.array([]), np.array([]), []
    results = []
    for res in process_file(fname, settings):
        results.append(res)
    cutoff = get_cutoff(results, FDR=1)
    print cutoff
    
    aa_mass = get_aa_mass(settings)

    j = 0
    for res in results:
        for e, (_, seq, note) in zip(res['e-values'], res['candidates']):
            e, seq, note, RT, spectrum = float(e), seq, note, utils.get_RT(res['spectrum']), res['spectrum']
            j += 1
            if note == 'd' or float(e) > cutoff:
                pass
            else:
                RTexp.append(RT)
                neutral_mass, _, _ = get_info(res['spectrum'], res, settings)
                massdif = np.append(massdif, (neutral_mass - mass.fast_mass(seq, charge=0, aa_mass=aa_mass)) / neutral_mass * 1e6)
                new_params = get_fragment_mass_tol(spectrum, seq, settings)
                fragstd = np.append(fragstd, new_params['fmt'])
                seqs.append(str(seq))

    print len(RTexp), 'top PSMs with 1% FDR'
    RC_def = achrom.RCs_yoshida
    xdict = {}
    for key, val in RC_def['aa'].items():
        xdict[key] = [val, None]
    RC_dict = achrom.get_RCs_vary_lcp(seqs, RTexp)
    RC_dict_new = dict()
    for key, val in RC_dict['aa'].items():
        xdict[key][1] = val
    a, b, _, _ = auxiliary.linear_regression([x[0] for x in xdict.values() if x[1] != None], [x[1] for x in xdict.values() if x[1] != None])
    for key, x in xdict.items():
        if x[1] == None:
            x[1] = x[0] * a + b
        RC_dict_new[key] = x[1]
    RC_dict['aa'] = RC_dict_new
    deltaRT = [rtexp - achrom.calculate_RT(pep, RC_dict, raise_no_mod=False)
            for rtexp, pep in zip([mean(x) for x in RTexp], seqs)]

    def condition(spectrum, cand, _):
        return 1 <= percentileofscore(deltaRT, utils.get_RT(spectrum)
                        - achrom.calculate_RT(cand, RC_dict)
                    ) <= 99

    step = FDbinSize(fragstd)
    lside, rside = 0, 1
    mt_h, _ = np.histogram(fragstd, bins=np.arange(lside, rside, step))

    for idx, mt in enumerate(mt_h):
        if mt == 0:
            mt_h = mt_h[:idx]
            break
    threshold = mt_h.size * step
    fragstd = fragstd[fragstd <= threshold]
    best_mt = max(fragstd[fragstd < scoreatpercentile(fragstd, 95)])
    
    best_par_mt = max(massdif[massdif < scoreatpercentile(massdif, 99.5)])

    print massdif
    print 'NEW FRAGMENT MASS TOLERANCE = %s' % (best_mt, )
    print 'NEW PARENT MASS TOLERANCE = %s' % (best_par_mt, )
    settings.set('search', 'product accuracy', best_mt)
    settings.set('search', 'precursor accuracy value', best_par_mt)
    settings.set('scoring', 'condition', condition)
    return settings
