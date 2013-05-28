from scipy.stats import percentileofscore
from pyteomics import achrom
from .main import *

def smart_filtering(fname, settings):
    conditions = settings.get('smart filtering', 'conditions')
    cutoff = settings.getfloat('smart filtering', 'initial e-value cut-off')
    # somehow derive the proper e-value cut-off
    # TODO
    settings = copy(settings)
    settings.set('misc', 'first stage', '')
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

def achrom_rt_filtering(fname, settings):
    cutoff = settings.getfloat('smart filtering', 'initial e-value cut-off')
    # somehow derive the proper e-value cut-off
    # TODO
    settings = copy(settings)
    settings.set('misc', 'first stage', '')
    RTexp, seqs = [], []
    for res in process_file(fname, settings):
        for e, (_, seq) in zip(res['e-values'], res['candidates']):
            if e > cutoff:
                break
            RTexp.append(utils.get_RT(res['spectrum']))
            seqs.append(seq)
    RC_dict = achrom.get_RCs_vary_lcp(seqs, RTexp)
    deltaRT = [rtexp - achrom.calculate_RT(pep, RC_dict, raise_no_mod=False)
            for rtexp, pep in zip(RTexp, seqs)]
    def condition(spectrum, cand, _):
        return 1 <= percentileofscore(deltaRT, utils.get_RT(spectrum)
                        - achrom.calculate_RT(cand, RC_dict)
                    ) <= 99
    settings.set('scoring', 'condition', condition)
    return settings
