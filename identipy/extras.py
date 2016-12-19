from scipy.stats import percentileofscore, scoreatpercentile
from pyteomics import achrom, auxiliary as aux, parser, mass
from main import *
from scoring import get_fragment_mass_tol
import numpy as np
from utils import get_info, get_aa_mass, get_output, get_enzyme, get_output
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


def get_subset(results, settings, fdr=0.01):
    """Filter results to given FDR using top 1 candidates"""
    out = get_output(results, settings)
    subset = aux.filter(out, key=lambda x: x['e-values'][0],
            is_decoy = lambda x: x['candidates'][0][2] == 'd',
            fdr=fdr)
    return subset

def optimization(fname, settings):
    settings = settings.copy()
    settings.set('misc', 'first stage', '')
    efc = settings.get('scoring', 'e-values for candidates')
    settings.set('scoring', 'e-values for candidates', 1)

    results = process_file(fname, settings)
    filtered = get_subset(results, settings, fdr=0.01)
    print len(filtered), 'PSMs with 1% FDR.'
    if len(filtered) < 50:
        if len(filtered) < 10:
            print 'OPTIMIZATION ABORTED'
            return settings
        else:
            functions = [precursor_mass_optimization, fragment_mass_optimization,
                    missed_cleavages_optimization]
    else:
        functions = [
                rt_filtering,
                precursor_mass_optimization, fragment_mass_optimization,
                charge_optimization, missed_cleavages_optimization]
    for func in functions:
        settings = func(filtered, settings)
    settings.set('scoring', 'e-values for candidates', efc)
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
    print 'NEW charges = %s:%s' % (mincharge, maxcharge)
    settings.set('search', 'maximum charge', maxcharge)
    settings.set('search', 'minimum charge', mincharge)
    return settings

def precursor_mass_optimization(results, settings):
    settings_nopime = settings.copy()
    settings_nopime.set('search', 'precursor isotope mass error', '0')
    settings_nopime.set('search', 'shifts', '0')
    results = get_output(results, settings_nopime)

    settings = settings.copy()
    massdif = np.array([res['candidates'][0][4]['mzdiff']['ppm'] for res in results])

    best_par_mt_l = min(massdif[massdif > scoreatpercentile(massdif, 0.5)])
    best_par_mt_r = max(massdif[massdif < scoreatpercentile(massdif, 99.5)])
    print 'NEW PARENT MASS TOLERANCE = %s:%s' % (best_par_mt_l, best_par_mt_r)
    settings.set('search', 'precursor accuracy left', -(best_par_mt_l))
    settings.set('search', 'precursor accuracy right', best_par_mt_r)
    settings.set('search', 'precursor accuracy unit', 'ppm')
    return settings

def missed_cleavages_optimization(results, settings):
    settings = settings.copy()
    missedcleavages = np.array([parser.num_sites(str(res['candidates'][0][1]), get_enzyme(str(settings.get('search', 'enzyme'))))
        for res in results])
    best_missedcleavages = missedcleavages.max()
    for mc in range(best_missedcleavages, -1, -1):
        if float(missedcleavages[missedcleavages > mc].size) / missedcleavages.size < 0.02:
            best_missedcleavages = mc
    print 'NEW miscleavages = %s' % (best_missedcleavages, )
    settings.set('search', 'number of missed cleavages', best_missedcleavages)
    return settings

def fragment_mass_optimization(results, settings):
    settings = settings.copy()
#    fragmassdif = np.array([get_fragment_mass_tol(res['spectrum'], str(res['candidates'][0][1]), settings)['fmt'] for res in results])
    fragmassdif = []
    for res in results:
        fragmassdif.extend(get_fragment_mass_tol(res['spectrum'], str(res['candidates'][0][1]), settings)['fmt'])
    fragmassdif = np.array(fragmassdif)
    step = FDbinSize(fragmassdif)
    lside, rside = 0, 1
    mt_h, _ = np.histogram(fragmassdif, bins=np.arange(lside, rside, step))

    for idx, mt in enumerate(mt_h):
        if mt == 0:
            mt_h = mt_h[:idx]
            break
#   threshold = mt_h.size * step
#   fragmassdif = fragmassdif[fragmassdif <= threshold]
    # best_frag_mt = max(fragmassdif[fragmassdif < scoreatpercentile(fragmassdif, 97.5)])
    best_frag_mt = scoreatpercentile(fragmassdif, 50) * 3

    print 'NEW FRAGMENT MASS TOLERANCE ppm = %s' % (best_frag_mt, )
    settings.set('search', 'product accuracy ppm', best_frag_mt)
    return settings


def rt_filtering(results, settings):
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
        print 'RT is missing. Skipping RT optimization.'
        return settings
    RC_def = achrom.RCs_gilar_rp
    xdict = {}
    for key, val in RC_def['aa'].items():
        xdict[key] = [val, None]
    RC_dict = achrom.get_RCs_vary_lcp(seqs, RTexp)
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
                    kk = ('-' + k[-2:]) if k[-1] == ']' else (k[-2:] + '-')
                else:
                    kk = k[:-1]
            elif len(k) > 1:
                kk = k
            print k, '->', kk
            RC_dict_new[v] = RC_dict_new[kk]

    RC_dict['aa'] = RC_dict_new

    print RC_dict
    rtexp = np.array([np.mean(x) for x in RTexp])
    rttheor = np.array([achrom.calculate_RT(pep, RC_dict, raise_no_mod=False)
        for pep in seqs])
    deltaRT = rtexp - rttheor
    print aux.linear_regression(rtexp, rttheor)
    print 'deltaRT percentiles:', scoreatpercentile(deltaRT, [1, 25, 50, 75, 99])
    
    h = FDbinSize(deltaRT)
    heights, edges = np.histogram(deltaRT, bins=np.arange(deltaRT.min(), deltaRT.max()+h, h))

    def condition(spectrum, cand, _):
        b = np.digitize(utils.get_RT(spectrum) - achrom.calculate_RT(list(cand), RC_dict),
                edges)

        return b and b < edges.size and heights[b-1]

    settings.set('scoring', 'condition', condition)
    return settings
