import re
from pyteomics import mass, electrochem as ec, auxiliary as aux, fasta, mzml, parser, mgf
import pandas as pd
from itertools import combinations, islice
from collections import defaultdict, Counter
import numpy as np
from multiprocessing import Queue, Process, cpu_count
import string
from copy import copy
try:
    from ConfigParser import RawConfigParser
except ImportError:
    from configparser import RawConfigParser
import tempfile
import os
import platform
import logging
import itertools as it
try:
    from lxml import etree
except ImportError:
    etree = None
from time import strftime
from os import path
logger = logging.getLogger(__name__)

try:
    from pyteomics import cmass
except ImportError:
    logger.warning('pyteomics.cythonize not found. It is highly recommended for good performance.')
    cmass = mass
try:
    from . import cparser
except ImportError:
    from . import customparser as cparser
from scipy.spatial import cKDTree
import pkg_resources
try:
    basestring
except NameError:
    basestring = (str, bytes)

default_tags = {
    'tmt10plex': {
        'tmt_126': 126.1277261,
        'tmt_127N': 127.1247610,
        'tmt_128C': 128.1344357,
        'tmt_129N': 129.1314706,
        'tmt_130C': 130.1411453,
        'tmt_131': 131.1381802,
        'tmt_127C': 127.1310809,
        'tmt_128N': 128.1281158,
        'tmt_129C': 129.1377905,
        'tmt_130N': 130.1348254
    },
    'tmt11plex': {
        'tmt_126': 126.1277261,
        'tmt_127N': 127.1247610,
        'tmt_127C': 127.1310809,
        'tmt_128N': 128.1281158,
        'tmt_128C': 128.1344357,
        'tmt_129N': 129.1314706,
        'tmt_129C': 129.1377905,
        'tmt_130N': 130.1348254,
        'tmt_130C': 130.1411453,
        'tmt_131': 131.1381802,
        'tmt_131C': 131.144999
    },
   'tmt_pro': {
        'tmt_126': 126.1277261,
        'tmt_127N': 127.1247610,
        'tmt_127C': 127.1310809,
        'tmt_128N': 128.1281158,
        'tmt_128C': 128.1344357,
        'tmt_129N': 129.1314706,
        'tmt_129C': 129.1377905,
        'tmt_130N': 130.1348254,
        'tmt_130C': 130.1411453,
        'tmt_131N': 131.1381802,
        'tmt_131C': 131.1445,
        'tmt_132N': 132.14153,
        'tmt_132C': 132.14785,
        'tmt_133N': 133.14489,
        'tmt_133C': 133.15121,
        'tmt_134N': 134.14824,
    },
    'tmt6plex': {
        'tmt_126': 126.1277261,
        'tmt_127N': 127.1247610,
        'tmt_128C': 128.1344357,
        'tmt_129N': 129.1314706,
        'tmt_130C': 130.1411453,
        'tmt_131': 131.1381802,
    }
}
default_tags['tmt16plex'] = default_tags['tmt_pro']


def get_tags(tags):
    logger.debug('Tags: %s', tags)
    if tags:
        if tags in default_tags:
            return default_tags[tags]
        else:
            ctags = dict()
            for tag in str(tags).split(','):
                for lbl, mss in tag.split(':'):
                    ctags[lbl] = float(mss)
            return ctags
    else:
        return tags


def get_child_for_mods(mods_str, settings, fixed=True, protein=False):
    if mods_str:
        for mod in re.split(r'[,;]\s*', mods_str):
            term = False
            if '-' not in mod:
                child_mod = etree.Element('aminoacid_modification')

                t = None
                if '[' in mod:
                    t = 'n'
                elif ']' in mod:
                    t = 'c'
                if t:
                    child_mod.set('protein_terminus' if protein else 'peptide_terminus', t)
                mod = mod.replace('[', '').replace(']', '')
                mod_label, mod_aa = parser._split_label(mod)
                mod_mass = mass.std_aa_mass.get(mod_aa, 0)
                mod_massdiff = settings.getfloat('modifications', mod_label)

                child_mod.set('aminoacid', mod_aa)
                child_mod.set('massdiff', str(mod_massdiff))
                child_mod.set('mass', str(mod_mass+mod_massdiff))
                child_mod.set('variable', 'Y' if not fixed else 'N')
                yield child_mod
            elif mod[0] == '-':
                term = 'c'
                mod_label = mod[1:]
                mod_term_mass = settings.getfloat('modifications', 'protein cterm cleavage')
            elif mod[-1] == '-':
                term = 'n'
                mod_label = mod[:-1]
                mod_term_mass = settings.getfloat('modifications', 'protein nterm cleavage')

            if term:
                mod_massdiff = settings.getfloat('modifications', mod_label)
                child_mod = etree.Element('terminal_modification')
                child_mod.set('terminus', term)
                child_mod.set('massdiff', str(mod_massdiff))
                child_mod.set('mass', str((mod_massdiff if not fixed else 0)+mod_term_mass))
                child_mod.set('variable', 'Y' if not fixed else 'N')
                yield child_mod


def custom_mass(sequence, nterm_mass, cterm_mass, **kwargs):
    return cmass.fast_mass(sequence, **kwargs) + (nterm_mass - 1.007825) + (cterm_mass - 17.002735)


def get_RCs(sequences, RTs, lcp=-0.21, term_aa=False, **kwargs):

    peptide_lengths = kwargs.get('lengths', np.log([len(peptide) for peptide in sequences]))
    peptide_dicts = sequences#[Counter(peptide) for peptide in sequences]

    detected_amino_acids = {aa for peptide_dict in peptide_dicts
                                for aa in peptide_dict}

    # Determine retention coefficients using multidimensional linear
    # regression.
    composition_array = []
    for idx, pdict in enumerate(peptide_dicts):
        loglen = peptide_lengths[idx]#np.log(parser.length(pdict))
        composition_array.append([pdict.get(aa, 0.) * (1. + lcp * loglen)
               for aa in detected_amino_acids] + [1.])

    # Add normalizing conditions for terminal retention coefficients. The
    # condition we are using here is quite arbitrary. It implies that the sum
    # of N- or C-terminal RCs minus the sum of corresponding internal RCs must
    # be equal to zero.
    if term_aa:
        for term_label in ['nterm', 'cterm']:
            normalizing_peptide = []
            for aa in detected_amino_acids:
                if aa.startswith(term_label):
                    normalizing_peptide.append(1.0)
                elif (term_label+aa) in detected_amino_acids:
                    normalizing_peptide.append(-1.0)
                else:
                    normalizing_peptide.append(0.0)
            normalizing_peptide.append(0.0)
            composition_array.append(normalizing_peptide)
            RTs.append(0.0)

    # Use least square linear regression.
    RCs, res, rank, s = np.linalg.lstsq(np.array(composition_array), np.array(RTs))

    # Remove normalizing elements from the RTs vector.
    if term_aa:
        for term_label in ['nterm', 'cterm']:
            RTs.pop()

    # Form output.
    RC_dict = {}
    RC_dict['aa'] = dict(
        zip(list(detected_amino_acids),
            RCs[:len(detected_amino_acids)]))
    RC_dict['aa'][parser.std_nterm] = 0.0
    RC_dict['aa'][parser.std_cterm] = 0.0
    RC_dict['const'] = RCs[len(detected_amino_acids)]
    RC_dict['lcp'] = lcp

    # Find remaining terminal RCs.
    if term_aa:
        for term_label in ['nterm', 'cterm']:
            # Check if there are terminal RCs remaining undefined.
            undefined_term_RCs = [aa for aa in RC_dict['aa']
                                if aa[1:5] != 'term' and term_label + aa not in RC_dict['aa']]
            if not undefined_term_RCs:
                continue

            # Find a linear relationship between internal and terminal RCs.
            defined_term_RCs = [aa for aa in RC_dict['aa']
                              if aa[1:5] != 'term' and term_label + aa in RC_dict['aa']]

            a, b, r, stderr = aux.linear_regression(
                [RC_dict['aa'][aa] for aa in defined_term_RCs],
                [RC_dict['aa'][term_label+aa] for aa in defined_term_RCs])

            # Define missing terminal RCs using this linear equation.
            for aa in undefined_term_RCs:
                RC_dict['aa'][term_label + aa] = a * RC_dict['aa'][aa] + b

    return RC_dict


def get_RCs_vary_lcp(sequences, RTs, term_aa=False, lcp_range=(-1.0, 1.0), **kwargs):

    labels = kwargs.get('labels')

    best_r = -1.1
    best_RC_dict = {}
    lcp_accuracy = kwargs.get('lcp_accuracy', 0.1)

    min_lcp = lcp_range[0]
    max_lcp = lcp_range[1]
    step = (max_lcp - min_lcp) / 10.0
    peptide_lengths = np.log([len(peptide) for peptide in sequences])
    peptide_dicts = [Counter(peptide) for peptide in sequences]
    while step > lcp_accuracy:
        lcp_grid = np.arange(min_lcp, max_lcp,
                                (max_lcp - min_lcp) / 10.0)
        for lcp in lcp_grid:
            RC_dict = get_RCs(peptide_dicts, RTs, lcp, term_aa, labels=labels, lengths=peptide_lengths)
            regression_coeffs = aux.linear_regression(
                RTs,
                [calculate_RT(peptide, RC_dict) for peptide in peptide_dicts])
            if regression_coeffs[2] > best_r:
                best_r = regression_coeffs[2]
                best_RC_dict = dict(RC_dict)
        min_lcp = best_RC_dict['lcp'] - step
        max_lcp = best_RC_dict['lcp'] + step
        step = (max_lcp - min_lcp) / 10.0

    return best_RC_dict


def calculate_RT(peptide, RC_dict, raise_no_mod=True):
    plen = len(peptide)
    peptide_dict = peptide
    RT = 0.0
    for aa in peptide_dict:
        if aa not in RC_dict['aa']:
            if len(aa) == 1:
                raise aux.PyteomicsError('No RC for residue "{}".'.format(aa))
            if (not raise_no_mod) and aa[-1] in RC_dict['aa']:
                RT += RC_dict['aa'][aa[-1]]
            else:
                raise aux.PyteomicsError(
                    'Residue "{0}" not found in RC_dict. '.format(aa) +
                    'Set raise_no_mod=False to ignore this error ' +
                    'and use the RC for "{0}"" instead.'.format(aa[-1]))
        else:
            RT += RC_dict['aa'][aa]

    length_correction_term = (
        1.0 + RC_dict.get('lcp', 0) * np.log(plen))
    RT *= length_correction_term

    RT += RC_dict.get('const', 0)

    return RT


_modchars = set(string.ascii_lowercase + string.digits)


def custom_split_label(mod):
    j = 0
    while mod[j] in _modchars:
        j += 1
    if j == 0:
        return mod[1:], '-', ']'
    if len(mod[j:]) > 1 and '[' in mod:
        return mod[:j], mod[j:].replace('[', ''), '['
    elif len(mod[j:]) > 1 and ']' in mod:
        return mod[:j], mod[j:].replace(']', ''), ']'
    elif len(mod[j:]) == 1:
        if mod.startswith('-'):
            return mod[:j], '-', ']'
        elif mod.endswith('-'):
            return mod[:j], '-', '['
        else:
            return mod[:j], mod[j:], ''


def iterate_spectra(fname):
    ftype = fname.rsplit('.', 1)[-1].lower()
    if ftype == 'mgf':
        with mgf.read(fname, read_charges=False, use_index=False) as f:
            for x in f:
                yield x
    elif ftype == 'mzml':
        with mzml.read(fname, use_index=False) as f:
            for x in f:
                if x['ms level'] > 1:
                    yield x
    else:
        raise ValueError('Unrecognized file type: {}'.format(ftype))


def iterate_and_preprocess(fname, params, settings):
    it = iterate_spectra(fname)
    n = settings.getint('performance', 'processes')
    return multimap(n, preprocess_spectrum, it, kwargs=params)


def is_decoy_function(settings):
    prefix = settings.get('input', 'decoy prefix').strip()
    infix = settings.get('input', 'decoy infix').strip()
    if infix:
        return lambda d: infix in d
    if prefix:
        return lambda d: d.startswith(prefix)
    logger.error('No decoy label specified. One of "decoy prefix" or "decoy infix" is needed.')


def peptide_gen(settings, clear_seen_peptides=False):
    if clear_seen_peptides:
        seen_target.clear()
        seen_decoy.clear()
    isdecoy = is_decoy_function(settings)
    enzyme = get_enzyme(settings.get('search', 'enzyme'))
    logger.debug('Using cleavage rule: %s', enzyme)
    semitryptic = settings.getint('search', 'semitryptic')
    mc = settings.getint('search', 'number of missed cleavages')
    minlen = settings.getint('search', 'peptide minimum length')
    maxlen = settings.getint('search', 'peptide maximum length')
    snp = settings.getint('search', 'snp')
    clip_M = settings.getboolean('search', 'clip N-terminal methionine')
    for prot in prot_gen(settings):
        for pep, pos in prot_peptides(prot[1], enzyme, mc, minlen, maxlen,
            is_decoy=isdecoy(prot[0]), snp=snp, desc=prot[0], semitryptic=semitryptic, position=True, clip_M=clip_M):
            term = ''
            if pos == 0:
                term += 'n'
            if pos + len(pep) == len(prot[1]):
                term += 'c'
            yield pep, term


def peptide_isoforms(settings, clear_seen_peptides=False):
    snp = settings.getint('search', 'snp')
    maxmods = settings.getint('modifications', 'maximum variable mods')
    leg = settings.get('misc', 'legend')
    pleg = settings.get('misc', 'plegend')
    logger.debug('leg: %s, pleg: %s', leg, pleg)
    punct = set(string.punctuation)
    nmods = [(p, mod[1], mod[2]) for p, mod in leg.items() if p in punct]
    pmods_n, pmods_c = [], []
    for p, mod in pleg.items():
        if p in punct:
            if mod[2] == '[':
                pmods_n.append((p, mod[1], mod[2]))
            if mod[2] == ']':
                pmods_c.append((p, mod[1], mod[2]))
    logger.debug('nmods: %s', nmods)
    logger.debug('pmods_n: %s', pmods_n)
    logger.debug('pmods_c: %s', pmods_c)
    aa_mass = get_aa_mass(settings)
    nterm_mass = settings.getfloat('modifications', 'protein nterm cleavage')
    cterm_mass = settings.getfloat('modifications', 'protein cterm cleavage')
    for peptide, term in peptide_gen(settings, clear_seen_peptides):
        mods = nmods[:]
        if 'n' in term:
            mods += pmods_n
        if 'c' in term:
            mods += pmods_c
        for form in (custom_isoforms(peptide, variable_mods=mods, maxmods=maxmods, snp=snp) if (nmods and maxmods) else [peptide, ]):
            if snp:
                if 'snp' not in form:
                    seqm = form
                    aachange_pos = False
                    snp_label = 'wild'
                else:
                    tmp = form.split('snp')
                    seqm = tmp[0] + tmp[1].split('at')[0].split('to')[-1] + tmp[2]
                    aachange_pos = len(tmp[0]) + 1
                    snp_label = tmp[1]
                aachange_pos = False
            else:
                seqm = form
                aachange_pos = False
                snp_label = False

            m = custom_mass(seqm, aa_mass=aa_mass, nterm_mass=nterm_mass, cterm_mass=cterm_mass)
            yield (seqm, aachange_pos, snp_label, m)


def prot_gen(settings):
    db = settings.get('input', 'database')
    # add_decoy = settings.getboolean('input', 'add decoy')
    # prefix = settings.get('input', 'decoy prefix')

    with fasta.read(db) as f:
        for p in f:
            yield p


def get_peptides(prot_seq, enzyme, mc, minlen, maxlen, semitryptic=False):
    peptides = cparser._cleave(prot_seq, enzyme, mc)
    for pep, startposition in peptides:
        plen = len(pep)
        if minlen <= plen <= maxlen:
            if not semitryptic:
                yield pep, startposition, plen
            else:
                for i in range(plen-minlen+1):
                    yield pep[i:], startposition + i, plen - i
                for i in range(1, plen-minlen+1, 1):
                    yield pep[:-i], startposition, plen - i


seen_target = set()
seen_decoy = set()
def prot_peptides(prot_seq, enzyme, mc, minlen, maxlen, is_decoy,
        dont_use_seen_peptides=False, snp=False, desc=False, position=False, semitryptic=False, clip_M=True):

    dont_use_fast_valid = parser.fast_valid(prot_seq)
    methionine_check = (clip_M and prot_seq[0] == 'M')
    if snp == 2:
        if desc:
            try:
                tmp = desc.split(' ')[0].split('|')
                pos = int(tmp[1]) - 1
                aach = tmp[2]
            except:
                desc = False
    # peptides = cparser._cleave(prot_seq, enzyme, mc)
    # for pep, startposition in peptides:
    #     plen = len(pep)
    for pep, startposition, plen in get_peptides(prot_seq, enzyme, mc, minlen, maxlen, semitryptic):
        loopcnt = 0
        if pep not in seen_target and pep not in seen_decoy and (dont_use_fast_valid or parser.fast_valid(pep)):
            loopcnt = 1
            if methionine_check and startposition == 0:
                if minlen <= plen - 2:
                    loopcnt = 3
                elif minlen <= plen - 1:
                    loopcnt = 2
        while loopcnt:
            f = pep[loopcnt-1:]
            if dont_use_seen_peptides:
                if snp == 1:
                    for ff, seq_new in custom_snp(f, startposition):
                        if not seq_new:
                            yield ff if not position else (ff, startposition)
                        else:
                            yield ff if not position else (ff, startposition)
                else:
                    yield f if not position else (f, startposition)
            else:
                if f not in seen_target and f not in seen_decoy:
                    if is_decoy:
                        seen_decoy.add(f)
                    else:
                        seen_target.add(f)
                    if snp == 1:
                        for ff, seq_new in custom_snp(f, startposition):
                            if not seq_new:
                                yield ff if not position else (ff, startposition)
                            if seq_new not in seen_decoy and seq_new not in seen_target:
                                yield ff if not position else (ff, startposition)
                    elif snp == 2:
                        if desc and startposition <= pos <= startposition + plen:
                            if len(aach) == 3 and aach[0] in parser.std_amino_acids and aach[2] in parser.std_amino_acids:
                                pos_diff = pos - startposition
                                f = f[:pos_diff] + 'snp%sto%sat%ssnp' % (aach.split('>')[0], aach.split('>')[-1], pos) + f[pos_diff+1:]
                                yield f if not position else (f, startposition)
                        else:
                            yield f if not position else (f, startposition)
                    else:
                        yield f if not position else (f, startposition)
            loopcnt -= 1


def custom_snp(peptide, startposition):
    yield peptide, None
    j = len(peptide) - 1
    while j >= 0:
        for aa in parser.std_amino_acids:
            if aa != 'L' and aa != peptide[j] and not (aa == 'I' and peptide[j] == 'L'):
                aa_label = 'snp%sto%sat%ssnp' % (peptide[j], aa, str(j + startposition))
                out = peptide[:j] + aa_label + peptide[j+1:], peptide[:j] + aa + peptide[j+1:]
                yield out
        j -= 1


def normalize_mods(sequence, settings):
    leg = settings.get('misc', 'legend')
    if leg:
        for char in string.punctuation:
            if char in leg:
                if leg[char][2] == ']' and leg[char][1] == '-':
                    sequence = sequence.replace(char, '-' + leg[char][0])
                else:
                    sequence = sequence.replace(char, ''.join(leg[char][:2]))
    return sequence


def custom_isoforms(peptide, variable_mods, maxmods=2, nterm=False, cterm=False, snp=False):
    if not variable_mods:
        yield peptide
    else:
        to_char = variable_mods[-1][0]
        from_char = variable_mods[-1][1]
        term = variable_mods[-1][2]
        sites = [s[0] for s in enumerate(peptide) if (not snp or (s[0] - 4 < 0 or peptide[s[0]-4:s[0]-1] != 'snp')) and (from_char == '-' or s[1] == from_char) and (not term or (term == '[' and s[0] == 0) or (term == ']' and s[0] == len(peptide)-1))]
        for m in range(maxmods+1):
            for comb in combinations(sites, m):
                flag = 0
                flag2 = 0
                tmpnterm = True if nterm else False
                tmpcterm = True if cterm else False
                v = ''
                cc_prev = 0
                for cc in comb:
                    if from_char == '-':
                        if term == '[' and not nterm:
                            flag2 = 1
                            v += to_char
                            tmpnterm = True
                        elif term == ']' and not cterm:
                            v = v + peptide[cc_prev:cc+1] + to_char
                            tmpcterm = True
                        else:
                            flag = 1
                    else:
                        v = v + peptide[cc_prev:cc] + to_char
                    if not flag2:
                        cc_prev = cc + 1
                if not flag:
                    v = v + peptide[cc_prev:]
                    for z in custom_isoforms(v, variable_mods[:-1], maxmods=maxmods - m, nterm=tmpnterm, cterm=tmpcterm, snp=snp):
                        yield z


def remove_precursor(mz_prec, spectrum, acc):
    mz = spectrum['m/z array']
    intens = spectrum['intensity array']
    idx = np.full(mz.size, True)
    i_l = mz.searchsorted(mz_prec - acc)#mz.size-2
    i_r = mz.searchsorted(mz_prec + acc, side='right')
    for i in range(i_l, i_r, 1):
        idx[i] = False
    spectrum['m/z array'] = mz[idx]
    spectrum['intensity array'] = intens[idx]


def deisotope(spectrum, acc, charge):
    #   acc = 0.3
    mz = spectrum['m/z array']
    intens = spectrum['intensity array']

    h = 1.0057
    i = mz.size-2
    skip = set()
    add = []
    while i >= 0:
        j = min(mz.size-1, mz.searchsorted(mz[i] + 1.5, side='right'))
        while j > i:
            if intens[i] > intens[j]:
                d = mz[j] - mz[i]
                if d > 1.5*h:
                    j -= 1
                    continue
                for z in range(1, charge+1):
                    if abs(d - 1./z) < acc:
                        skip.add(j)
                        if z > 1:
    #                         skip.add(i)
                            add.append((i, z))
            j -= 1
        i -= 1
    ix = np.delete(np.arange(mz.size, dtype=int), list(skip))
    newmz, newint = [], []
    for i, z in add:
        newmz.append(mz[i]*z - (z-1)*h)
        newint.append(intens[i])
    #   print len(skip), len(add)
    mz = np.hstack((mz[ix], newmz))
    intens = np.hstack((intens[ix], newint))
    spectrum['m/z array'] = mz
    spectrum['intensity array'] = intens


def preprocess_spectrum(spectrum, kwargs):
    spectrum = copy(spectrum)
    maxpeaks = kwargs['maxpeaks']
    minpeaks = kwargs['minpeaks']
    dynrange = kwargs['dynrange']
    acc = kwargs['acc']
    tags = kwargs['tags']

    if 'm/z array' not in spectrum:
        return None

    _, states = get_expmass(spectrum, kwargs)
    if not states:
        return None

    if tags:
        # TODO optimize performance
        max_mass_label_val = max(tags.values()) + 1.0
        tmp_idx = np.nonzero(spectrum['m/z array'] <= max_mass_label_val)
        tags_res = defaultdict(float)
        for tmt_label, tmt_mass in tags.items():
            for t_m, t_i in zip(spectrum['m/z array'][tmp_idx], spectrum['intensity array'][tmp_idx]):
                if abs(t_m - tmt_mass) / tmt_mass <= 1e-5:
                    tags_res[tmt_label] += t_i
        for tmt_label, tmt_intensity in tags_res.items():
            spectrum[tmt_label] = tmt_intensity

    if kwargs['deisotope']:
        dacc = kwargs['dacc']
        deisotope(spectrum, dacc, states[-1])

    mz_prec, _ = get_expmass(spectrum, kwargs)
    remove_precursor(mz_prec, spectrum, acc)

    mz = spectrum['m/z array']

    idx = np.nonzero(mz >= kwargs['min_mz'])
    spectrum['intensity array'] = spectrum['intensity array'][idx]
    mz = mz[idx]
    spectrum['intensity array'] = spectrum['intensity array'].astype(np.float32)

    if minpeaks and spectrum['intensity array'].size < minpeaks:
        return None

    spectrum['intensity array'] = spectrum['intensity array'].astype(np.float32)

    if dynrange:
        i = spectrum['intensity array'] > spectrum['intensity array'].max(
                ) / dynrange
        spectrum['intensity array'] = spectrum['intensity array'][i]
        mz = mz[i]

    if maxpeaks and minpeaks > maxpeaks:
        raise ValueError('minpeaks > maxpeaks: {} and {}'.format(
            minpeaks, maxpeaks))
    if maxpeaks and spectrum['intensity array'].size > maxpeaks:
        i = np.argsort(spectrum['intensity array'])[-maxpeaks:]
        j = np.argsort(mz[i])
        spectrum['intensity array'] = spectrum['intensity array'][i][j]
        mz = mz[i][j]

    spectrum['m/z array'] = mz

    if minpeaks and spectrum['intensity array'].size < minpeaks:
        return None

    spectrum['Isum'] = spectrum['intensity array'].sum()

    tmp2 = dict()
    tmp = spectrum['m/z array'] / acc
    tmp = tmp.astype(int)
    for idx, mt in enumerate(tmp):
        i_val = spectrum['intensity array'][idx] / spectrum['Isum']
        for mz_val_int in (mt-1, mt, mt+1):
            if mz_val_int not in tmp2:
                tmp2[mz_val_int] = i_val
            else:
                tmp2[mz_val_int] = max(i_val, tmp2[mz_val_int])
    tmp = np.concatenate((tmp, tmp-1, tmp+1))
    spectrum['fastset'] = set(tmp.tolist())
    spectrum['RT'] = get_RT(spectrum)
    spectrum['idict'] = tmp2

    spectrum['__KDTree'] = cKDTree(spectrum['m/z array'].reshape((spectrum['m/z array'].size, 1)))

    return spectrum


def relative(unit):
    if unit == 'ppm':
        return True
    elif unit in {'Th', 'Da', 'amu'}:
        return False
    else:
        raise ValueError('Unrecognized precursor accuracy unit: ' + unit)


def set_mod_dict(settings):
    mods = settings.get('modifications', 'variable')
    pmods = settings.get('modifications', 'protein variable')

    settings.set('modifications', 'variable_original', mods)
    settings.set('modifications', 'protein_original', pmods)
    i = None
    if isinstance(mods, basestring):
        mods = mods.strip(' ,')
        mod_dict = {}
        legend = {}

        if mods:
            mods = [custom_split_label(l) for l in re.split(r',\s*', mods)]
            mods.sort(key=lambda x: len(x[0]), reverse=True)
            for i, (mod, char) in enumerate(zip(mods, string.punctuation), 1):
                legend[''.join(mod)] = char
                legend[char] = mod
            assert all(len(m) == 3 for m in mods), 'unmodified residue given'
            for mod, aa, term in mods:
                mod_dict.setdefault(mod, []).append(aa)
        settings.set('modifications', 'variable', mod_dict)
        logger.info('Setting legend: %s', legend)
        settings.set('misc', 'legend', legend)

    if isinstance(pmods, basestring):
        plegend = {}
        pmod_dict = {}
        if pmods:
            pmods = [custom_split_label(l) for l in re.split(r',\s*', pmods)]
            pmods.sort(key=lambda x: len(x[0]), reverse=True)
            for mod, char in zip(pmods, string.punctuation[i:]):
                plegend[''.join(mod)] = char
                plegend[char] = mod
            assert all(len(m) == 3 for m in pmods), 'unmodified residue given'
            for mod, aa, term in pmods:
                pmod_dict.setdefault(mod, []).append(aa)
        settings.set('modifications', 'protein variable', pmod_dict)
        mod_dict.update(pmod_dict)
        settings.set('modifications', 'variable', mod_dict)
        settings.set('misc', 'plegend', plegend)
        logger.info('Setting plegend: %s', plegend)


def get_enzyme(enzyme):
    if enzyme in parser.expasy_rules:
        return parser.expasy_rules[enzyme]
    else:
        try:
            enzyme = convert_tandem_cleave_rule_to_regexp(enzyme)
            return enzyme
        except Exception as e:
            logger.debug('Exception parsing cleavage rule %s: %s', enzyme, e.args[0])
            return enzyme


def convert_tandem_cleave_rule_to_regexp(cleavage_rule):

    def get_sense(c_term_rule, n_term_rule):
        if '{' in c_term_rule:
            return 'N'
        elif '{' in n_term_rule:
            return 'C'
        else:
            if len(c_term_rule) <= len(n_term_rule):
                return 'C'
            else:
                return 'N'

    def get_cut(cut, no_cut):
        aminoacids = set(parser.std_amino_acids)
        cut = ''.join(aminoacids & set(cut))
        if '{' in no_cut:
            no_cut = ''.join(aminoacids & set(no_cut))
            return cut, no_cut
        else:
            no_cut = ''.join(set(parser.std_amino_acids) - set(no_cut))
            return cut, no_cut

    out_rules = []
    for protease in cleavage_rule.split(','):
        protease = protease.replace('X', ''.join(parser.std_amino_acids))
        c_term_rule, n_term_rule = protease.split('|')
        sense = get_sense(c_term_rule, n_term_rule)
        if sense == 'C':
            cut, no_cut = get_cut(c_term_rule, n_term_rule)
        else:
            cut, no_cut = get_cut(n_term_rule, c_term_rule)

        if no_cut:
            if sense == 'C':
                out_rules.append('([%s](?=[^%s]))' % (cut, no_cut))
            else:
                out_rules.append('([^%s](?=[%s]))' % (no_cut, cut))
        else:
            if sense == 'C':
                out_rules.append('([%s])' % (cut, ))
            else:
                out_rules.append('(?=[%s])' % (cut, ))
    return '|'.join(out_rules)


class CustomRawConfigParser(RawConfigParser, object):
    def get(self, section, option, **kwargs):
        val = super(CustomRawConfigParser, self).get(section, option)
        if isinstance(val, basestring):
            if section == 'search' and option == 'enzyme':
                return val.split('|class')[0]
            return val[::-1].split('|', 1)[-1][::-1]
        return val

    def get_choices(self, section, option):
        val = super(CustomRawConfigParser, self).get(section, option)
        if isinstance(val, basestring) and len(val.split('|')) > 1:
            return val[::-1].split('|', 1)[0][::-1]
        else:
            return ''

    def copy(self):
        new_config = CustomRawConfigParser()
        for section in self.sections():
            new_config.add_section(section)
            for name, value in self.items(section):
                new_config.set(section, name, value)
        return new_config


def find_nearest(array, value):
    return (np.abs(np.array(array) - value)).argmin()


def _charge_params(settings):
    params = {}
    params['maxcharge'] = settings.getint('search', 'maximum charge') or None
    params['mincharge'] = settings.getint('search', 'minimum charge') or None
    if settings.has_option('search', 'minimum unknown charge') and settings.getint('search', 'minimum unknown charge'):
        params['min_ucharge'] = max(settings.getint('search', 'minimum unknown charge'), params['mincharge'])
    else:
        params['min_ucharge'] = params['mincharge']
    if settings.has_option('search', 'maximum unknown charge') and settings.getint('search', 'maximum unknown charge'):
        params['max_ucharge'] = min(settings.getint('search', 'maximum unknown charge'), params['maxcharge'])
    else:
        params['max_ucharge'] = params['maxcharge']
    return params


def get_info(spectrum, result, settings, aa_mass=None):
    'Returns neutral mass, charge state and retention time of the top candidate'
    if not aa_mass:
        aa_mass = get_aa_mass(settings)
    RT = spectrum['RT']#get_RT(spectrum)

    params = _charge_params(settings)

    masses, states = zip(*neutral_masses(spectrum, params))
    # idx = find_nearest(masses, cmass.fast_mass(str(result['candidates'][0][1]), aa_mass=aa_mass))


    nterm_mass = settings.getfloat('modifications', 'protein nterm cleavage')
    cterm_mass = settings.getfloat('modifications', 'protein cterm cleavage')

    idx = find_nearest(masses, custom_mass(str(result['candidates'][0][1]), aa_mass=aa_mass, nterm_mass=nterm_mass, cterm_mass=cterm_mass))
    return (masses[idx], states[idx], RT)


def reshape_theor_spectrum(peaks):
    for k in peaks.keys():
        marr = np.array(peaks[k])
        n = marr.size
        peaks[k] = marr.reshape((n, 1))
    return peaks


ion_shift_dict = {
    'a': 46.00547930326002,
    'b': 18.010564683699954,
    'c': 0.984015582689949,
    'x': -25.979264555419945,
    'y': 0.0,
    'z': 17.026549101010005,
}


def calc_ions_from_neutral_mass(peptide, nm, ion_type, charge, aa_mass, cterm_mass, nterm_mass):
    if ion_type in 'abc':
        nmi = nm - aa_mass[peptide[-1]] - ion_shift_dict[ion_type] - (cterm_mass - 17.002735)
    else:
        nmi = nm - aa_mass[peptide[0]] - ion_shift_dict[ion_type] - (nterm_mass - 1.007825)
    return (nmi + 1.0072764667700085 * charge) / charge


def check_n_term(ion_type):
    return (ion_type[0] == 'b' or ion_type[0] == 'a' or ion_type[0] == 'c')


def get_n_ions(peptide, maxmass, pl, charge, k_aa_mass):
    tmp = [maxmass, ]
    for i in range(1, pl):
        tmp.append(tmp[-1] - k_aa_mass[peptide[-i-1]]/charge)
    return tmp


def get_c_ions(peptide, maxmass, pl, charge, k_aa_mass):
    tmp = [maxmass, ]
    for i in range(pl-2, -1, -1):
        tmp.append(tmp[-1] - k_aa_mass[peptide[-(i+2)]]/charge)
    return tmp


def theor_spectrum(peptide, acc_frag, nterm_mass, cterm_mass, types=('b', 'y'), maxcharge=None, reshape=False, **kwargs):
    peaks = {}
    theoretical_set = dict()
    if 'nm' in kwargs:
        nm = kwargs['nm']
    else:
        nm = custom_mass(peptide, aa_mass=kwargs['aa_mass'], nterm_mass = nterm_mass, cterm_mass = cterm_mass)
    pl = len(peptide) - 1
    if not maxcharge:
        maxcharge = 1 + int(ec.charge(peptide, pH=2))
    for charge in range(1, maxcharge + 1):
        for ion_type in types:
            nterminal = check_n_term(ion_type)
            if nterminal:
                maxmass = calc_ions_from_neutral_mass(peptide, nm, ion_type=ion_type, charge=charge,
                                aa_mass=kwargs['aa_mass'], cterm_mass=cterm_mass, nterm_mass=nterm_mass)
                marr = get_n_ions(peptide, maxmass, pl, charge, kwargs['aa_mass'])
            else:
                maxmass = calc_ions_from_neutral_mass(peptide, nm, ion_type=ion_type, charge=charge,
                                aa_mass=kwargs['aa_mass'], cterm_mass=cterm_mass, nterm_mass=nterm_mass)
                marr = get_c_ions(peptide, maxmass, pl, charge, kwargs['aa_mass'])

            tmp = [int(x / acc_frag) for x in marr]
            if ion_type in theoretical_set:
                theoretical_set[ion_type].extend(tmp)
            else:
                theoretical_set[ion_type] = tmp

            if reshape:
                marr = np.array(marr)
                n = marr.size
                marr = marr.reshape((n, 1))
            peaks[ion_type, charge] = marr
    return peaks, theoretical_set


def get_expmass(spectrum, kwargs):
    maxcharge = kwargs['maxcharge'] or None
    mincharge = kwargs['mincharge'] or None
    min_ucharge = kwargs['min_ucharge']
    max_ucharge = kwargs['max_ucharge']

    if 'params' in spectrum:
        exp_mass = spectrum['params']['pepmass'][0]
        charge = spectrum['params'].get('charge')
    else:
        ion = spectrum['precursorList']['precursor'][
                0]['selectedIonList']['selectedIon'][0]
        charge = ion.get('charge state')
        if charge is not None: charge = [int(charge)]
        exp_mass = ion['selected ion m/z']

    if isinstance(charge, str):
        states = [s for s in aux._parse_charge(charge, True)
                if (mincharge is None or s >= mincharge) and (maxcharge is None or s <= maxcharge)]
    elif charge is None:
        states = list(range(min_ucharge, 1 + max_ucharge))
    else:
        states = [c for c in charge if
            (mincharge is None or c >= mincharge) and (maxcharge is None or c <= maxcharge)]
    states.sort()
    return exp_mass, states


def neutral_masses(spectrum, params):
    exp_mass, states = get_expmass(spectrum, params)
    return zip((c * (exp_mass - mass.nist_mass['H+'][0][0])
            for c in states), states)


@aux.memoize(10)
def import_(name):
    """Import a function by name: module.function or
    module.submodule.function, etc. By default trying to find
    function name in identipy.scoring module.
    Return the function object."""

    try:
        mod, f = name.rsplit('.', 1)
        return getattr(__import__(mod, fromlist=[f]), f)
    except Exception as e:
        logger.error('%s', e)
        return getattr(__import__('identipy.scoring', fromlist=[name]), name)


def get_aa_mass(settings):
    if settings.has_option('misc', 'aa_mass'):
        return settings.get('misc', 'aa_mass')
    aa_mass = mass.std_aa_mass.copy()
    aa_mass['-'] = 0.0
    for k, v in settings.items('modifications'):
        if k not in {'fixed', 'variable', 'variable_original', 'protein variable', 'protein_original'}:
            aa_mass[k] = float(v)
    fmods = settings.get('modifications', 'fixed')
    if fmods:
        for mod in re.split(r'[,;]\s*', fmods):
            if '-' not in mod:
                m, aa = parser._split_label(mod)
                aa_mass[aa] += settings.getfloat('modifications', m)
    vmods = settings.get('modifications', 'variable')
    if vmods:
        leg = settings.get('misc', 'legend')
        for p in string.punctuation:
            if p in leg:
                mod, aa, term = leg[p]
                if term == ']' and aa == '-':
                    aa_mass[p] = aa_mass[mod] + aa_mass[aa]
                    aa_mass[aa+mod] = aa_mass[mod] + aa_mass[aa]
                else:
                    aa_mass[p] = aa_mass[mod] + aa_mass[aa]
                    aa_mass[mod+aa] = aa_mass[mod] + aa_mass[aa]
    pmods = settings.get('modifications', 'protein variable')
    if pmods:
        leg = settings.get('misc', 'plegend')
        for p in string.punctuation:
            if p in leg:
                mod, aa, term = leg[p]
                if term == ']' and aa == '-':
                    aa_mass[p] = aa_mass[mod] + aa_mass[aa]
                    aa_mass[aa+mod] = aa_mass[mod] + aa_mass[aa]
                else:
                    aa_mass[p] = aa_mass[mod] + aa_mass[aa]
                    aa_mass[mod+aa] = aa_mass[mod] + aa_mass[aa]
    return aa_mass


def multimap(n, func, it, global_data, best_res_in=False, best_res_raw_in=False, best_peptides=False, **kw):
    global best_res


    rel = kw['rel']
    nterm_mass = kw.get('nterm_mass')
    cterm_mass = kw.get('cterm_mass')
    acc_l = kw['acc_l']
    acc_r = kw['acc_r']

    shifts_and_pime = kw['sapime']

    if best_res_in:
        best_res = deepcopy(best_res_in)
        best_res_raw = deepcopy(best_res_raw_in)
    else:
        best_res = {}
        best_res_raw = {}
    best_res_pep = {}

    if n == 0:
        try:
            n = cpu_count()
        except NotImplementedError:
            n = 1
    if n == 1 or platform.system() == 'Windows':
        cnt1 = 0
        for s in it:
            cnt1 += 1
            if cnt1 % 10000 == 0:
                logger.debug(cnt1)
            result = func(s, best_res, global_data[0], **kw)
            if result:
                for x in result:
                    peptide, m, snp_label, res = x

                    for score, spec_t, c, info in res:
                        if -score <= best_res.get(spec_t, 0):
                            best_res_raw[spec_t] = [peptide, m, snp_label, score, spec_t, c, info]
                            best_res[spec_t] = -score
        return best_res_raw, best_res

    else:

        def worker(qout, start, end, global_data_local):
            # maxval = len(qin)
            # start = 0

            new_best_res = {}
            new_best_res_raw = {}
            best_pep_res = {}

            while start < end:
                item = qin[start]
                result = func(item, best_res, global_data_local, **kw)

                if result:
                    for x in result:
                        peptide, m, snp_label, res = x

                        for score, spec_t, c, info in res:
                            if -score <= new_best_res.get(spec_t, best_res.get(spec_t, 0)):
                                new_best_res[spec_t] = -score
                                best_res[spec_t] = -score
                                new_best_res_raw[spec_t] = [peptide, m, snp_label, score, spec_t, c, info]
                start += 1
            qout.put(new_best_res_raw)
            qout.put(None)
        qsize = kw.pop('qsize')
        qout = Queue(qsize)
        count = 0

        global qin

        while True:
            qint = list(islice(it, 5000000))
            if not len(qint):
                break

            qin = []
            for seqm, aachange_pos, snp_label, m in qint:
                qin.append((seqm, aachange_pos, snp_label, m))
            qin = sorted(qin, key=lambda x: x[3])
            qin_masses = np.array([z[3] for z in qin])

            procs = []
            for proc_num in range(n):

                min_mass = min(global_data[proc_num]['nmasses'])
                max_mass = max(global_data[proc_num]['nmasses'])
                if rel:
                    dm_l = acc_l * max_mass / 1.0e6
                    dm_r = acc_r * max_mass / 1.0e6
                elif not rel:
                    dm_l = acc_l
                    dm_r = acc_r
                dm_l -= min(shifts_and_pime)
                dm_r += max(shifts_and_pime)
                start = qin_masses.searchsorted(min_mass + dm_l)
                end = qin_masses.searchsorted(max_mass + dm_r)

                p = Process(target=worker, args=(qout, start, end, global_data[proc_num]))
                p.start()
                procs.append(p)

            count = len(qin)

            for _ in range(n):
                logger.debug('%s %s', _, len(best_res_pep))
                for item in iter(qout.get, None):
                    for k, v in item.items():
                        if -v[3] <= best_res.get(k, 0):
                            best_res_raw[k] = v
                            best_res[k] = -v[3]
                            best_res_pep[k] = v[0]
                logger.debug('%s %s', _, len(best_res_pep))

            logger.debug('HERE1')

            for p in procs:
                p.join()

            logger.debug('HERE2')

        logger.info(len(best_res_pep))
        return best_res_raw, best_res


def allow_all(*args):
    return True


def get_RT(spectrum):
    """Return scan retention time in seconds"""
    # MGF
    if 'params' in spectrum:
        try:
            return float(spectrum['params']['rtinseconds'])
        except:
            try:
                return float(spectrum['params']['title'].split(',')[-1].strip().split()[0])
            except:
                try:
                    return 60 * np.average([float(x) for x in spectrum['params']['title'].split('lution from: ')[-1].split(' period:')[0].split(' to ')])
                except:
                    return 0
    # mzML
    try:
        rt = spectrum['scanList']['scan'][0]['scan start time']
        try:
            if rt.unit_info == 'second':
                return float(rt)
            else:
                return float(rt * 60)
        except AttributeError:
            return float(rt)
    except KeyError:
        return 0


def get_title(spectrum):
    if 'params' in spectrum:
        return spectrum['params']['title']
    else:
        return spectrum['id']


def get_precursor_mz(spectrum):
    try:
        return spectrum['params']['pepmass'][0]
    except:
        return spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']


def is_db_target_only(settings):
    db = settings.get('input', 'database')
    isdecoy = is_decoy_function(settings)
    balance = 0
    for prot in fasta.read(db):
        if isdecoy(prot[0]):
            balance -= 1
        else:
            balance += 1
    return bool(balance)


def get_shifts_and_pime(settings):
    # try:
    #     fast_first_stage = settings.getint('misc', 'fast first stage')
    # except:
    #     fast_first_stage = 0
    pime = settings.getint('search', 'precursor isotope mass error')
    # if fast_first_stage:
    #     pime = 0
    shifts =[float(x) for x in settings.get('search', 'shifts').split(',')]
    dM = mass.nist_mass['C'][13][0] - mass.nist_mass['C'][12][0]
    shifts_and_pime = shifts[:]
    for i in range(pime):
        shifts_and_pime += [x + (i + 1) * dM for x in shifts]
    return shifts_and_pime


def build_pept_prot(settings, results):
    mc = settings.getint('search', 'number of missed cleavages')
    minlen = settings.getint('search', 'peptide minimum length')
    maxlen = settings.getint('search', 'peptide maximum length')
    isdecoy = is_decoy_function(settings)
    clip_M = settings.getboolean('search', 'clip N-terminal methionine')

    snp = settings.getint('search', 'snp')
    pept_prot = {}
    prots = {}
    peptides = set()
    pept_neighbors = {}
    pept_ntts = {}
    enzyme = settings.get('search', 'enzyme')
    semitryptic = settings.getint('search', 'semitryptic')
    for x in results:
        peptides.update(re.sub(r'[^A-Z]', '', normalize_mods(x['candidates'][i][1], settings)) for i in range(
            1 or len(x['candidates'])))
    seen_target.clear()
    seen_decoy.clear()
    enzyme_rule = get_enzyme(enzyme)
    for desc, prot in prot_gen(settings):
        dbinfo = desc.split(' ')[0]
        prots[dbinfo] = desc
        if semitryptic:
            cl_positions = set(z for z in it.chain([x.end() for x in re.finditer(enzyme_rule, prot)],
                   [0, 1, len(prot)]))
        for pep, startposition in prot_peptides(prot, enzyme_rule, mc, minlen, maxlen, isdecoy(desc),
            dont_use_seen_peptides=True, snp=snp, desc=desc, position=True, semitryptic=semitryptic, clip_M=clip_M):
            if snp:
                if 'snp' not in pep:
                    seqm = pep
                else:
                    tmp = pep.split('snp')
                    seqm = tmp[0] + tmp[1].split('at')[0].split('to')[-1] + tmp[2]
            else:
                seqm = pep

            if seqm in peptides:
                pept_ntts.setdefault(seqm, {})
                pept_neighbors.setdefault(seqm, {})
                pept_neighbors[seqm][dbinfo] = (prot[startposition - 1] if startposition != 0 else '-',
                        prot[startposition + len(seqm)] if startposition + len(seqm) < len(prot) else '-')

                if not semitryptic:
                    pept_prot.setdefault(seqm, []).append(dbinfo)
                    pept_ntts[seqm][dbinfo] = 2
                else:
                    ntt = (startposition in cl_positions) + ((startposition + len(seqm)) in cl_positions)
                    pept_ntts[seqm][dbinfo] = ntt
                    pept_prot.setdefault(seqm, []).append(dbinfo)

    return pept_prot, prots, pept_neighbors, pept_ntts


def get_outpath(inputfile, settings, suffix):
    outpath = settings.get('output', 'path')
    filename = os.path.join(outpath, os.path.splitext(os.path.basename(inputfile))[0] + os.path.extsep + suffix)
    return filename


def write_pepxml(inputfile, settings, results):
    outpath = settings.get('output', 'path')
    logger.debug('Output path: %s', outpath)

    set_mod_dict(settings)

    enzyme = settings.get('search', 'enzyme')
    search_engine = 'IdentiPy'
    database = settings.get('input', 'database')
    missed_cleavages = settings.getint('search', 'number of missed cleavages')
    fmods = settings.get('modifications', 'fixed')
    snp = settings.getint('search', 'snp')
    nterm_mass = settings.getfloat('modifications', 'protein nterm cleavage')
    cterm_mass = settings.getfloat('modifications', 'protein cterm cleavage')
    tags = get_tags(settings.get('output', 'tags'))

    nterm_fixed = 0
    cterm_fixed = 0

    for mod in re.split(r'[,;]\s*', fmods):
        if mod.startswith('-'):
            cterm_fixed = settings.getfloat('modifications', 'protein cterm cleavage')
        elif mod.endswith('-'):
            nterm_fixed = settings.getfloat('modifications', 'protein nterm cleavage')

    filename = get_outpath(inputfile, settings, 'pep.xml')
    with open(filename, 'wb') as output:
        logger.info('Writing %s ...', filename)
        line1 = b'<?xml version="1.0" encoding="UTF-8"?>\n\
        <?xml-stylesheet type="text/xsl" href="pepXML_std.xsl"?>\n'
        output.write(line1)

        base_name, ftype = path.splitext(inputfile)
        ftype = ftype.lower()

        root = etree.Element('msms_pipeline_analysis')
        root.set("date", strftime("%Y:%m:%d:%H:%M:%S"))
        root.set("summary_xml", '')
        root.set("xmlns", 'http://regis-web.systemsbiology.net/pepXML')
        # TODO
        #root.set("xmlns:xsi", 'http://www.w3.org/2001/XMLSchema-instance')
        #root.set("xsi:schemaLocation", 'http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v117.xsd')

        child1 = etree.Element('msms_run_summary')
        child1.set("base_name", base_name)
        child1.set("search_engine", search_engine)
        child1.set("raw_data_type", "raw")  # ?

        if ftype == '.mgf':
            child1.set("raw_data", ".mgf")
        elif ftype == '.mzml':
            child1.set("raw_data", ".mzML")
        else:
            child1.set("raw_data", ".?")
        root.append(child1)

        child2 = etree.Element('sample_enzyme')
        child2.set('name', enzyme)
        child1.append(child2)

        child3 = etree.Element('specificity')
        child3.set("cut", "KR")
        child3.set("no_cut", "P")
        child3.set("sense", "C")

        child2.append(child3)

        child4 = etree.Element('search_summary')
        child4.set('base_name', base_name)
        child4.set('search_engine', search_engine)
        child4.set("search_engine_version", get_version())
        child4.set('precursor_mass_type', 'monoisotopic')
        child4.set('fragment_mass_type', 'monoisotopic')
        child4.set('search_id', '1')

        for child_mod in get_child_for_mods(settings.get('modifications', 'fixed'), settings, fixed=True):
            child4.append(child_mod)
        for child_mod in get_child_for_mods(settings.get('modifications', 'variable_original'), settings, fixed=False):
            child4.append(child_mod)
        for child_mod in get_child_for_mods(settings.get('modifications', 'protein_original'), settings, fixed=False, protein=True):
            child4.append(child_mod)

        child1.append(child4)

        child5 = etree.Element('search_database')
        child5.set('local_path', database)
        child5.set('type', 'AA')

        child4.append(copy(child5))

        child5 = etree.Element('enzymatic_search_constraint')
        child5.set('enzyme', enzyme)
        child5.set('max_num_internal_cleavages', str(missed_cleavages))
        child5.set('min_number_termini', '2')

        child4.append(copy(child5))

        results = [x for x in results if x['candidates'].size]
#       results = list(get_output(results, settings))
        logger.info('Accumulated results: %s', len(results))
        pept_prot, prots, pept_neighbors, pept_ntts = build_pept_prot(settings, results)
        if settings.has_option('misc', 'aa_mass'):
            aa_mass = settings.get('misc', 'aa_mass')
        else:
            aa_mass = get_aa_mass(settings)
        vmods = set()
        variablemods = settings.get('modifications', 'variable')
        if variablemods:
            for k, v in variablemods.items():
                for aa in v:
                    vmods.add(k + aa)
                    vmods.add(aa + k)

        leg = {}
        if settings.has_option('misc', 'legend'):
            leg = settings.get('misc', 'legend')
        if settings.has_option('misc', 'plegend'):
            leg.update(settings.get('misc', 'plegend'))

        ntermcleavage = settings.getfloat('modifications', 'protein nterm cleavage')
        ctermcleavage = settings.getfloat('modifications', 'protein cterm cleavage')

        for idx, result in enumerate(results):
            if result['candidates'].size:
                tmp = etree.Element('spectrum_query')
                spectrum = result['spectrum']
                tmp.set('spectrum', get_title(spectrum))
                tmp.set('spectrumNativeID', get_title(spectrum))
                tmp.set('start_scan', str(idx))  # ???
                tmp.set('end_scan', str(idx))  # ???
                tmp.set('index', str(idx))  # ???

                neutral_mass, charge_state, RT = get_info(spectrum, result, settings, aa_mass)
                tmp.set('precursor_neutral_mass', str(neutral_mass))
                tmp.set('assumed_charge', str(int(charge_state)))
                if RT:
                    tmp.set('retention_time_sec', str(RT))

                tmp2 = etree.Element('search_result')
                result['candidates'] = result['candidates'][:len(result['e-values'])]

                flag = 1
                for i, candidate in enumerate(result['candidates']):
                    match = candidate[4]['match']
                    if match is None:
                        break
                    tmp3 = etree.Element('search_hit')
                    tmp3.set('hit_rank', str(i + 1))
                    mod_sequence = str(candidate[1])
                    mod_sequence = normalize_mods(mod_sequence, settings)
                    sequence = re.sub(r'[^A-Z]', '', mod_sequence)
                    if sequence not in pept_prot:
                        flag = 0
                        logger.error('Unaccounted sequence! %s (%s)', sequence, mod_sequence)
                        break
                    else:
                        tmp3.set('peptide', sequence)

                        proteins = pept_prot[re.sub(r'[^A-Z]', '', sequence)]

                        tmp3.set('protein', prots[proteins[0]].split(' ', 1)[0] + (('_' + candidate[7]) if snp else ''))
                        try:
                            protein_descr = prots[proteins[0]].split(' ', 1)[1]
                        except:
                            protein_descr = ''

                        neighbors = pept_neighbors.get(sequence, {}).get(proteins[0], ('-', '-'))

                        tmp3.set('peptide_prev_aa', neighbors[0])
                        tmp3.set('peptide_next_aa', neighbors[1])
                        tmp3.set('protein_descr', protein_descr)

                        num_tot_proteins = len(proteins)
                        tmp3.set('num_tot_proteins', str(num_tot_proteins))
                        tmp3.set('num_matched_ions', str(sum(v.sum() for v in match.values())))
                        tmp3.set('tot_num_ions', str((len(sequence) - 1) * 2))
                        neutral_mass_theor = custom_mass(sequence, aa_mass=aa_mass, nterm_mass=nterm_mass, cterm_mass=cterm_mass)
                        # neutral_mass_theor = cmass.fast_mass(sequence, aa_mass=aa_mass)
                        tmp3.set('calc_neutral_pep_mass', str(neutral_mass_theor))
                        tmp3.set('massdiff', str(candidate[4]['mzdiff']['Da']))
                        tmp3.set('num_tol_term', str(pept_ntts.get(sequence, {}).get(proteins[0], '?')))
                        tmp3.set('num_missed_cleavages', str(len(parser.cleave(sequence, get_enzyme(enzyme), 0)) - 1))
                        tmp3.set('is_rejected', '0')  # ???

                        if num_tot_proteins > 1 and (not snp or 'wild' not in prots[proteins[0]].split(' ', 1)[0]):
                            for prot in proteins[1:]:
                                tmp4 = etree.Element('alternative_protein')
                                tmp4.set('protein', prots[prot].split(' ', 1)[0] + (('_' + candidate[7]) if snp else ''))
                                try:
                                    protein_descr = prots[prot].split(' ', 1)[1]
                                except:
                                    protein_descr = ''
                                tmp4.set('protein_descr', protein_descr)
                                neighbors = pept_neighbors.get(sequence, {}).get(prot, ('-', '-'))
                                tmp4.set('peptide_prev_aa', neighbors[0])
                                tmp4.set('peptide_next_aa', neighbors[1])
                                tmp4.set('num_tol_term', str(pept_ntts.get(sequence, {}).get(prot, '?')))
                                tmp3.append(copy(tmp4))

                        labels = parser.std_labels + [la[:-1] if la[-1] == '[' else '-' + la[:-2] if la[-1] == ']' else la for la in leg if len(la) > 1]
#                       logger.debug('Known labels: %s', labels)
                        try:
                            aalist = parser.parse(mod_sequence, labels=labels)
                        except Exception as e:
                            logger.debug('Problematic sequence: %s\n%s', mod_sequence, e)
                            aalist = [a[::-1] for a in parser.parse(mod_sequence[::-1], labels=labels)][::-1]
                        tmp4 = etree.Element('modification_info')
                        ntermmod = 0

                        if nterm_fixed:
                            tmp4.set('mod_nterm_mass', str(nterm_fixed))
                        if cterm_fixed:
                            tmp4.set('mod_cterm_mass', str(cterm_fixed))

                        for idx, aminoacid in enumerate(aalist):
                            if aminoacid in fmods or aminoacid in vmods:
                                if aminoacid.endswith('-') and idx == 0:
                                    ntermmod = 1
                                    tmp4.set('mod_nterm_mass', str(str(aa_mass.get(aminoacid) + ntermcleavage)))
                                elif aminoacid.startswith('-') and idx == len(aalist) - 1:
                                    tmp4.set('mod_cterm_mass', str(aa_mass.get(aminoacid) + ctermcleavage))
                                else:
                                    tmp5 = etree.Element('mod_aminoacid_mass')
                                    tmp5.set('position', str(idx + 1 - ntermmod))
                                    tmp5.set('mass', str(aa_mass.get(aminoacid)))
                                    tmp4.append(copy(tmp5))
                        tmp3.append(copy(tmp4))

                        if 'RNHS' in candidate[4]:
                            tmp4 = etree.Element('search_score')
                            tmp4.set('name', 'hyperscore')
                            tmp4.set('value', str(candidate[4]['RNHS']))
                            tmp3.append(copy(tmp4))

                            tmp4 = etree.Element('search_score')
                            tmp4.set('name', 'expect')
                            tmp4.set('value', str(1./candidate[4]['RNHS']))
                            tmp3.append(copy(tmp4))

                        else:
                            tmp4 = etree.Element('search_score')
                            tmp4.set('name', 'hyperscore')
                            tmp4.set('value', str(candidate[0]))
                            tmp3.append(copy(tmp4))

                            tmp4 = etree.Element('search_score')
                            tmp4.set('name', 'expect')
                            tmp4.set('value', str(result['e-values'][i]))
                            tmp3.append(copy(tmp4))

                        tmp4 = etree.Element('search_score')
                        tmp4.set('name', 'sumI')
                        tmp4.set('value', str(candidate[5]))
                        tmp3.append(copy(tmp4))

                        tmp4 = etree.Element('search_score')
                        tmp4.set('name', 'fragmentMT')
                        tmp4.set('value', str(candidate[6]))
                        tmp3.append(copy(tmp4))

                        tmp4 = etree.Element('search_score')
                        tmp4.set('name', 'nextscore_std')
                        tmp4.set('value', str(candidate[8]))
                        tmp3.append(copy(tmp4))

                        if 'params' in spectrum:
                            if 'isowidthdiff' in spectrum['params']:
                                tmp4 = etree.Element('search_score')
                                tmp4.set('name', 'ISOWIDTHDIFF')
                                tmp4.set('value', str(spectrum['params'].get('isowidthdiff', 0)))
                                tmp3.append(copy(tmp4))

                            if 'rtwidth' in spectrum['params']:
                                tmp4 = etree.Element('search_score')
                                tmp4.set('name', 'RTwidth')
                                tmp4.set('value', str(spectrum['params'].get('rtwidth', 0)))
                                tmp3.append(copy(tmp4))

                            if 'ms1intensity' in spectrum['params']:
                                tmp4 = etree.Element('search_score')
                                tmp4.set('name', 'MS1Intensity')
                                tmp4.set('value', str(spectrum['params'].get('ms1intensity', 0)))
                                tmp3.append(copy(tmp4))

                            if 'pif' in spectrum['params']:
                                tmp4 = etree.Element('search_score')
                                tmp4.set('name', 'PIF')
                                tmp4.set('value', str(spectrum['params'].get('pif', -3)))
                                tmp3.append(copy(tmp4))

                            if 'sulfur' in spectrum['params']:
                                tmp4 = etree.Element('search_score')
                                tmp4.set('name', 'sulfur')
                                tmp4.set('value', str(spectrum['params'].get('sulfur', -1)))
                                tmp3.append(copy(tmp4))

                            if 'ionmobility' in spectrum['params']:
                                tmp4 = etree.Element('search_score')
                                tmp4.set('name', 'ionmobility')
                                tmp4.set('value', str(spectrum['params'].get('ionmobility', 0)))
                                tmp3.append(copy(tmp4))

                        if tags:
                            for tag_label in tags.keys():
                                tmp4 = etree.Element('search_score')
                                tmp4.set('name', 'tag_' + tag_label)
                                tmp4.set('value', str(spectrum.get(tag_label, 0)))
                                tmp3.append(copy(tmp4))


                        for k, v in match.items():
                            tmp4 = etree.Element('search_score')
                            tmp4.set('name', 'matched_{}{}_ions'.format(*k))
                            tmp4.set('value', str(v.sum()))
                            tmp3.append(copy(tmp4))

                        tmp2.append(copy(tmp3))
                if flag:
                    tmp.append(copy(tmp2))
                    child1.append(copy(tmp))

        s = etree.tostring(root, pretty_print=True)
        output.write(s)


def write_csv(inputfile, settings, results):
    df = dataframe(inputfile, settings, results)
    if df is None:
        logger.info('No results to write. File not created.')
        return

    sep = settings.get('output', 'separator')
    of = settings.get('output', 'format').lower()
    if not sep:
        sep = ',' if of == 'csv' else '\t'
    fname = get_outpath(inputfile, settings, of)
    logger.info('Writing %s ...', fname)
    df.to_csv(fname, index=False, sep=sep)


def dataframe(inputfile, settings, results):
    #   results = list(get_output(results, settings))
    results = list(results)
    if not results:
        return None

    logger.info('Accumulated results: %s', len(results))
    #   ensure_decoy(settings)
    set_mod_dict(settings)
    fmods = settings.get('modifications', 'fixed')
    pept_prot, prots, pept_neighbors, pept_ntts = build_pept_prot(settings, results)
    if settings.has_option('misc', 'aa_mass'):
        aa_mass = settings.get('misc', 'aa_mass')
    else:
        aa_mass = get_aa_mass(settings)

    nterm_mass = settings.getfloat('modifications', 'protein nterm cleavage')
    cterm_mass = settings.getfloat('modifications', 'protein cterm cleavage')

    vmods = set()
    variablemods = settings.get('modifications', 'variable')
    if variablemods:
        for k, v in variablemods.items():
            for aa in v:
                vmods.add(k + aa)
                vmods.add(aa + k)

    leg = {}
    if settings.has_option('misc', 'legend'):
        leg = settings.get('misc', 'legend')

    enzyme = settings.get('search', 'enzyme')
    snp = settings.getint('search', 'snp')
    columns = ['Title', 'Assumed charge', 'RT', 'Rank', 'Matched ions', 'Total ions', 'Calculated mass',
                'Mass difference', 'Missed cleavages', 'Proteins', '# proteins', 'Sequence', 'Modified sequence',
                'Hyperscore', 'Expect', 'sumI', 'fragmentMT']
    rows = []
    for result in results:
        if result['candidates'].size:
            row = []
            spectrum = result['spectrum']
            row.append(get_title(spectrum))
            neutral_mass, charge_state, RT = get_info(spectrum, result, settings, aa_mass)
            row.append(charge_state)
            row.append(RT)
            result['candidates'] = result['candidates'][:len(result['e-values'])]

            flag = 1
            for i, candidate in enumerate(result['candidates'], 1):
                match = candidate[4]['match']
                if match is None: break
                row.append(i)
                mod_sequence = normalize_mods(candidate[1], settings)

                sequence = re.sub(r'[^A-Z]', '', mod_sequence)
                if sequence not in pept_prot:
                    flag = 0
                    logger.error('Unaccounted sequence! %s (%s)', sequence, mod_sequence)
                    break
                else:
                    allproteins = pept_prot[re.sub(r'[^A-Z]', '', sequence)]

                    row.append(sum(v.sum() for v in match.values()))
                    row.append((len(sequence) - 1) * 2)
                    neutral_mass_theor = custom_mass(sequence, aa_mass=aa_mass, nterm_mass = nterm_mass, cterm_mass = cterm_mass)
                    row.append(neutral_mass_theor)
                    row.append(candidate[4]['mzdiff']['Da'])
                    row.append(parser.num_sites(sequence, get_enzyme(enzyme)))

                    proteins = [allproteins[0]]
                    if len(allproteins) > 1:
                        if snp:
                            wilds = any('wild' in prots[p].split(' ', 1)[0] for p in allproteins)
                        for prot in allproteins[1:]:
                            d = prots[prot].split(' ', 1)[0]
                            if (not snp or not wilds or 'wild' in d):
                                proteins.append(prot)

                    row.append(';'.join(proteins))
                    row.append(len(proteins))

                    row.append(sequence)
                    if fmods:
                        for mod in re.split(r'[,;]\s*', fmods):
                            if '-' not in mod:
                                m, aa = parser._split_label(mod)
                                mod_sequence = mod_sequence.replace(aa, m+aa)
                            elif mod[0] == '-':
                                mod_sequence = mod_sequence + mod
                            elif mod[-1] == '-':
                                mod_sequence = mod + mod_sequence
                    row.append(mod_sequence)

                    row.append(candidate[0])
                    row.append(result['e-values'][i-1])
                    row.append(candidate[5])
                    row.append(candidate[6])

                    rows.append(row)
    df = pd.DataFrame(rows)
    df.columns = columns
    return df


def write_pickle(inputfile, settings, results):
    results = list(results)
    logger.info('Accumulated results: %s', len(results))
    try:
        import cPickle as pickle
    except ImportError:
        import pickle
    filename = get_outpath(inputfile, settings, 'pickle')
    with open(filename, 'wb') as output:
        pickle.dump((inputfile, settings, results), output, -1)


def write_output(inputfile, settings, results):
    formats = {'pepxml': write_pepxml, 'csv': write_csv, 'tsv': write_csv, 'pickle': write_pickle}
    of = settings.get('output', 'format')
    writer = formats[re.sub(r'[^a-z]', '', of.lower())]

    if settings.has_option('output', 'path'):
        outd = settings.get('output', 'path')
        if not os.path.isdir(outd):
            logger.info('Creating %s ...', outd)
            os.makedirs(outd)
    else:
        outpath = os.path.dirname(inputfile)
        settings.set('output', 'path', outpath)

    return writer(inputfile, settings, results)


def demix_chimeric(path_to_features, path_to_mzml, demixing=False, calc_PIF=True):

    basename_mzml = os.path.splitext(path.basename(path_to_mzml))[0]

    if path_to_features:
        df1 = pd.read_table(path_to_features)
        df1 = df1.rename(columns=lambda x: x.strip())
        df1 = df1[df1['nIsotopes'] >= 2]
        logger.info(df1.shape)
    else:
        df1 = None

    mzs = []
    RTs = []
    ionmobs = []
    chs = []
    titles = []
    ms2_map = {}
    isolation_window_left = False
    isolation_window_right = False

    cur_ms1 = False
    mass_acc = 20
    for a in mzml.read(path_to_mzml):
        if a['ms level'] == 1:
            cur_ms1 = a
        elif a['ms level'] == 2:
            logger.debug('PROCESSING RT: ' + str(a['scanList']['scan'][0]['scan start time']) )
        if a['ms level'] == 2 and 'm/z array' in a.keys() and len(a['m/z array']) > 0:
            # if :
            if not isolation_window_left:
                isolation_window_left = float(a['precursorList']['precursor'][0]['isolationWindow']['isolation window lower offset'])
                isolation_window_right = float(a['precursorList']['precursor'][0]['isolationWindow']['isolation window upper offset'])
            pepmass = float(a['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
            RT = float(a['scanList']['scan'][0]['scan start time'])
            # logger.info(a['mean inverse reduced ion mobility array'])
            try:
                ion_mob = float(a['mean inverse reduced ion mobility array'][0])
            except:
                # logger.info('missing ion mob')
                ion_mob = 0.0
            try:
                ch = int(a['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
            except:
                ch = 0

            if calc_PIF:
                if not cur_ms1:
                    calc_PIF = False
                    # logger.info('Missing MS1 spectra in mzML, turning off PIF calculation')
                elif not isolation_window_left:
                    calc_PIF = False
                    # logger.info('Missing isolation window info in mzML, turning off PIF calculation')
                else:
                    intensity_full_ms2 = 0
                    intensity_precursor = 0

                    idx_l = cur_ms1['m/z array'].searchsorted(pepmass - isolation_window_left)
                    idx_r = cur_ms1['m/z array'].searchsorted(pepmass + isolation_window_right)

                    if not ch:
                        tch = 2
                    else:
                        tch = ch

                    abs_error = pepmass * mass_acc * 1e-6
                    for mz, intensity in zip(cur_ms1['m/z array'][idx_l:idx_r], cur_ms1['intensity array'][idx_l:idx_r]):
                        if any(abs(mz - (pepmass + (k * 1.007825) / tch)) <= abs_error for k in [-2, -1, 0, 1, 2, 3, 4]):
                            intensity_precursor += intensity
                        intensity_full_ms2 += intensity
                    if intensity_full_ms2:
                        PIF = intensity_precursor / intensity_full_ms2 * 100
                    else:
                        PIF = -1
                    a['PIF'] = PIF


            title = a['id']
            mzs.append(pepmass)
            RTs.append(RT)
            ionmobs.append(ion_mob)
            chs.append(ch)
            titles.append(title)
            ms2_map[title] = a


    mzs = np.array(mzs)
    RTs = np.array(RTs)
    ionmobs = np.array(ionmobs)
    chs = np.array(chs)
    titles = np.array(titles)
    idx = np.argsort(mzs)
    mzs = mzs[idx]
    RTs = RTs[idx]
    ionmobs = ionmobs[idx]
    chs = chs[idx]
    titles = titles[idx]

    if not df1 is None:
        if 'ion_mobility' not in df1.columns:
            df1['ion_mobility'] = 0
        if 'sulfur' not in df1.columns:
            df1['sulfur'] = 0
        df1['MSMS'] = df1.apply(findMSMS, axis=1, args = (isolation_window_left, isolation_window_right, mzs, RTs, titles, ionmobs))
        df1['MSMS_accurate'] = df1.apply(findMSMS_accurate, axis=1, args = (mzs, RTs, titles, ionmobs, chs))
        # print(df1['MSMS_accurate'])

    outmgf_name = os.path.splitext(path_to_mzml)[0] + '_identipy' + os.extsep + 'mgf'
    outmgf = open(outmgf_name, 'w')

    t_i = 1
    f_i = 0

    added_MSMS = set()

    if demixing:

        for z in df1[['mz', 'rtApex', 'charge', 'intensityApex', 'MSMS', 'MSMS_accurate', 'rtStart', 'rtEnd', 'ion_mobility', 'sulfur']].values:
            mz, RT, ch, Intensity, ttls, ttl_ac, rt_ll, rt_rr, ion_mob, sulfur = z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], z[9]
            if ttls:
                f_i += 1
                for ttl in ttls:
                    if ttl in ttl_ac:
                        added_MSMS.add(ttl)
                    mz_arr, I_arr = ms2_map[ttl]['m/z array'], ms2_map[ttl]['intensity array']
                    PIF = ms2_map[ttl].get('PIF', -2)
                    pepmass = float(ms2_map[ttl]['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
                    t_i_orig = ms2_map[ttl]['index']
                    outmgf.write('BEGIN IONS\n')
                    outmgf.write('TITLE=%s.%d.%d.%d.%d\n' % (basename_mzml, t_i_orig, t_i, t_i, ch))
                    outmgf.write('RTINSECONDS=%f\n' % (RT * 60, ))
                    outmgf.write('PEPMASS=%f %f\n' % (mz, Intensity))
                    outmgf.write('CHARGE=%d+\n' % (ch, ))
                    outmgf.write('ISOWIDTHDIFF=%f\n' % (mz - pepmass, ))
                    outmgf.write('RTwidth=%f\n' % (rt_rr - rt_ll, ))
                    outmgf.write('MS1Intensity=%f\n' % (Intensity, ))
                    outmgf.write('PIF=%f\n' % (PIF, ))
                    outmgf.write('IonMobility=%f\n' % (ion_mob, ))
                    outmgf.write('Sulfur=%f\n' % (sulfur, ))
                    for mz_val, I_val in zip(mz_arr, I_arr):
                        outmgf.write('%f %f\n' % (mz_val, I_val))
                    outmgf.write('END IONS\n\n')
                    t_i += 1

        for k in ms2_map:
            if k not in added_MSMS:
                f_i += 1
                a = ms2_map[k]
                mz_arr, I_arr = a['m/z array'], a['intensity array']
                PIF = a.get('PIF', -2)
                mz = float(a['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
                RT = float(a['scanList']['scan'][0]['scan start time'])
                t_i_orig = a['index']
                try:
                    ch = int(a['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
                except:
                    ch = ''
                outmgf.write('BEGIN IONS\n')
                outmgf.write('TITLE=%s.%d.%d.%d.%s\n' % (basename_mzml, t_i_orig, t_i, t_i, str(ch)))
                outmgf.write('RTINSECONDS=%f\n' % (RT * 60, ))
                outmgf.write('PEPMASS=%f %f\n' % (mz, 0))
                if ch:
                    outmgf.write('CHARGE=%d+\n' % (ch, ))
                outmgf.write('ISOWIDTHDIFF=%f\n' % (0.0, ))
                outmgf.write('RTwidth=%f\n' % (0.0, ))
                outmgf.write('MS1Intensity=%f\n' % (0.0, ))
                outmgf.write('PIF=%f\n' % (PIF, ))
                outmgf.write('IonMobility=%f\n' % (0.0, ))
                outmgf.write('Sulfur=%f\n' % (-1.0, ))
                for mz_val, I_val in zip(mz_arr, I_arr):
                    outmgf.write('%f %f\n' % (mz_val, I_val))
                outmgf.write('END IONS\n\n')
                t_i += 1

    else:

        MS2_acc_map = {}

        if not df1 is None:
            for z in df1[['mz', 'rtApex', 'charge', 'intensityApex', 'MSMS', 'MSMS_accurate', 'rtStart', 'rtEnd', 'ion_mobility', 'sulfur']].values:
                ttl_ac = z[5]
                for ttl in ttl_ac:
                    if ttl not in MS2_acc_map:
                        MS2_acc_map[ttl] = z
                    else:
                        if MS2_acc_map[ttl][3] < z[3]:
                            MS2_acc_map[ttl] = z

        # print(MS2_acc_map)

        for k in ms2_map:
            a = ms2_map[k]

            if k in MS2_acc_map:
                # print('HERE, ok')
                z = MS2_acc_map[k]
                mz, RT, ch, Intensity, ttls, ttl_ac, rt_ll, rt_rr, ion_mob, sulfur = z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], z[9]
            else:
                mz = float(a['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
                RT = float(a['scanList']['scan'][0]['scan start time'])
                try:
                    ch = int(a['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
                except:
                    ch = 0
                Intensity, ttls, ttl_ac, rt_ll, rt_rr, ion_mob, sulfur = 0, 0, 0, 0, 0, 0, 0


            mz_arr, I_arr = a['m/z array'], a['intensity array']
            PIF = a.get('PIF', -2)
#             mz = float(a['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
#             RT = float(a['scanList']['scan'][0]['scan start time'])
            t_i_orig = a['index']
            outmgf.write('BEGIN IONS\n')
            outmgf.write('TITLE=%s.%d.%d.%d.%s\n' % (basename_mzml, t_i_orig, t_i, t_i, str(ch)))
            outmgf.write('RTINSECONDS=%f\n' % (RT * 60, ))
            outmgf.write('PEPMASS=%f %f\n' % (mz, Intensity))
            if ch:
                outmgf.write('CHARGE=%d+\n' % (ch, ))
            outmgf.write('ISOWIDTHDIFF=%f\n' % 0)
            outmgf.write('RTwidth=%f\n' % (rt_rr - rt_ll, ))
            outmgf.write('MS1Intensity=%f\n' % (Intensity, ))
            outmgf.write('PIF=%f\n' % (PIF, ))
            outmgf.write('IonMobility=%f\n' % (ion_mob, ))
            outmgf.write('Sulfur=%f\n' % (sulfur, ))
            for mz_val, I_val in zip(mz_arr, I_arr):
                outmgf.write('%f %f\n' % (mz_val, I_val))
            outmgf.write('END IONS\n\n')
            t_i += 1
    outmgf.close()

    return outmgf_name


def findMSMS(raw, isolation_window_left, isolation_window_right, mzs, RTs, titles, ionmobs):
    out = []
    isotope_fix = raw['nIsotopes'] / raw['charge']
    mz = raw['mz']
    RT_l = raw['rtStart']
    RT_r = raw['rtEnd']
    ion_mob_p = raw['ion_mobility']
    # There is no error below: -right and +left!
    id_l = mzs.searchsorted(mz - isolation_window_right)
    id_r = mzs.searchsorted(mz + isolation_window_left + isotope_fix, side='right')
    for idx, RT in enumerate(RTs[id_l:id_r]):
        if RT_l <= RT <= RT_r:
            if abs(ionmobs[id_l+idx] - ion_mob_p) <= 0.1:
                out.append(titles[id_l+idx])
    if len(out):
        return out
    else:
        return None


def findMSMS_accurate(raw, mzs, RTs, titles, ionmobs, chs):
    out = set()
    acc = 10
    mz = raw['mz']
    ch = raw['charge']
    RT_l = raw['rtStart']
    RT_r = raw['rtEnd']
    ion_mob_p = raw['ion_mobility']
    acc_rel = mz * acc * 1e-6
    id_l = mzs.searchsorted(mz - acc_rel)
    id_r = mzs.searchsorted(mz + acc_rel, side='right')
    for idx, RT in enumerate(RTs[id_l:id_r]):
        if RT_l <= RT <= RT_r:
            if abs(ionmobs[id_l+idx] - ion_mob_p) <= 0.1:
                ch_msms = chs[id_l+idx]
                if not ch_msms or ch_msms == ch:
                    out.add(titles[id_l+idx])
            # return True
    # return False
    return out


def generate_database(settings, outname=None):
    add_decoy = settings.getboolean('input', 'add decoy')
    prefix = settings.get('input', 'decoy prefix')
    infix = settings.get('input', 'decoy infix')
    if infix and add_decoy:
        if not prefix:
            prefix = infix
        logger.warning('infix is specified with "add decoy" = True. Generated decoys will have PREFIX %s', prefix)
    mode = settings.get('input', 'decoy method')
    db = settings.get('input', 'database')
    target_only = is_db_target_only(settings)
    if add_decoy and target_only:
        gdbname = outname or settings.get('output', 'generated database')
        if gdbname:
            ft = open(gdbname, 'w')
        else:
            ft = tempfile.NamedTemporaryFile(mode='w', delete=False)
        fasta.write_decoy_db(db, ft, mode=mode, prefix=prefix)
        ft.flush()
        settings.set('input', 'database', ft.name)
        settings.set('input', 'add decoy', 'no')
        logger.debug('Generated database: %s (isfile = %s)', ft.name, os.path.isfile(ft.name))
        return ft.name
    else:
        logger.debug('Skipping database generation. add_decoy = %s, target_only = %s', add_decoy, target_only)


def get_version():
    return pkg_resources.get_distribution('identipy').version
