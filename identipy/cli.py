from . import main, utils
import argparse
import string

def get_label(modmass, labels):
    abt = string.ascii_lowercase
    abt_l = len(abt) - 1
    if modmass in labels:
        return labels[modmass], labels, 0
    else:
        labels[modmass] = abt[labels['i']] + abt[labels['j']] + abt[labels['k']]
        if labels['k'] > abt_l:
            labels['k'] = 0
            labels['j'] += 1
        if labels['j'] > abt_l:
            labels['j'] = 0
            labels['i'] += 1
        return labels[modmass], labels, 1

def _update(settings, section, name, value):
    if value is not None:
        settings.set(section, name, value)

def run():
    parser = argparse.ArgumentParser(
        description='Search proteins using LC-MS/MS spectra',
        epilog='''

    Example usage
    -------------
    $ identipy input.mgf -db human.fasta
    -------------
    ''',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('file',     help='input .mzML or .mgf file with MS/MS spectra')
    parser.add_argument('-db',      help='path to protein fasta file')
    parser.add_argument('-cfg',     help='path to file with parameters')
    parser.add_argument('-punit',   help='precursor mass tolerance unit. Can be ppm or Da', default='ppm', type=str)
    parser.add_argument('-ptol',    help='precursor mass tolerance', default=10.0, type=float)
    parser.add_argument('-lptol',   help='*left precursor mass tolerance', default=False, type=float)
    parser.add_argument('-rptol',   help='*right precursor mass tolerance', default=False, type=float)
    parser.add_argument('-funit',   help='fragment mass tolerance unit. Can be ppm or Da', default='Da')
    parser.add_argument('-ftol',    help='fragment mass tolerance', default=0.3, type=float)
    parser.add_argument('-fminmz',  help='fragment min m/z', default=100.0, type=float)
    parser.add_argument('-lmin',    help='min length of peptides', default=5, type=int)
    parser.add_argument('-lmax',    help='max length of peptides', default=55, type=int)
    parser.add_argument('-massmin', help='min mass of peptides', default=300, type=int)
    parser.add_argument('-massmax', help='max mass of peptides', default=4000, type=int)
    parser.add_argument('-e',       help='cleavage rule in quotes!. X!Tandem style for cleavage rules', default='[RK]|{P}')
    parser.add_argument('-mc',      help='number of missed cleavages', default=0, type=int)
    parser.add_argument('-cmin',    help='min precursor charge', default=1, type=int)
    parser.add_argument('-cmax',    help='max precursor charge', default=9, type=int)
    parser.add_argument('-cumin',   help='min unknown precursor charge', default=0, type=int)
    parser.add_argument('-cumax',   help='max unknown precursor charge', default=0, type=int)
    parser.add_argument('-ime',     help='precursor isotope mass error. The parent ion\
     mass tolerance is expanded by opening up multiple tolerance windows centered\
      on the given number of 13C isotope peaks for a peptide.', default=0, type=int)
    parser.add_argument('-shifts',  help='shifts. example: 0,16.000,23.000,12', default='0')
    parser.add_argument('-snp',     help='1 means make SNP changes for ALL peptides', default=0, type=int)
    parser.add_argument('-mm',      help='number of minimum matched ions', default=4, type=int)
    parser.add_argument('-ad',      help='add decoy', default=0, type=int)
    parser.add_argument('-prefix',  help='decoy prefix', default='DECOY_')
    parser.add_argument('-method',  help='reverse or random', default='reverse')
    parser.add_argument('-deis',    help='use MS/MS deisotoping. yes or no', default='no')
    parser.add_argument('-deistol', help='deisotope mass accuracy', default=0.3, type=float)
    parser.add_argument('-score',   help='used score. Can be RNHS, hyperscore or morpheusscore', default='RNHS')
    parser.add_argument('-minp',    help='minumum peaks in MS/MS spectra', default=4, type=int)
    parser.add_argument('-maxp',    help='maximum peaks in MS/MS spectra', default=100, type=int)
    parser.add_argument('-dyn',     help='dynamic range', default=100.0, type=float)
    parser.add_argument('-mfc',     help='maximum fragment charge', default=1, type=int)
    parser.add_argument('-nproc',   help='number of processes. 0 means auto', default=0, type=int)
    parser.add_argument('-maxmods', help='maximum variable mods per sequence', default=2, type=int)
    parser.add_argument('-ncleave', help='protein nterm cleavage', default=1.007825, type=float)
    parser.add_argument('-ccleave', help='protein cterm cleavage', default=17.002735, type=float)
    parser.add_argument('-fmods',   help='fixed modifications. in mass1@aminoacid1,mass2@aminoacid2 format', default='57.021464@C')
    parser.add_argument('-vmods',   help='variable modifications. in mass1@aminoacid1,mass2@aminoacid2 format', default='')

    args = vars(parser.parse_args())
    if args['cfg']:
        settings = main.settings(args['cfg'])
    else:
        settings = main.settings()

    labels = {'i': 0, 'j': 0, 'k': 0}


    if args['fmods']:
        for mod in args['fmods'].split(','):
            modmass, modaa = mod.split('@')
            lbl, labels, flag = get_label(modmass, labels)
            settings.set('modifications', 'fixed', lbl + modaa)
            if flag:
                settings.set('modifications', lbl, modmass)

    if args['vmods']:
        for mod in args['vmods'].split(','):
            modmass, modaa = mod.split('@')
            lbl, labels, flag = get_label(modmass, labels)
            settings.set('modifications', 'variable', lbl + modaa)
            if flag:
                settings.set('modifications', lbl, modmass)

    _update(settings, 'input', 'database', args['db'])
    _update(settings, 'search', 'precursor accuracy unit', args['punit'])
    _update(settings, 'search', 'precursor accuracy left', (args['ptol'] if not args['lptol'] else args['lptol']))
    _update(settings, 'search', 'precursor accuracy right', (args['ptol'] if not args['rptol'] else args['rptol']))
    _update(settings, 'search', 'product accuracy unit', args['funit'])
    _update(settings, 'search', 'product accuracy', args['ftol'])
    _update(settings, 'search', 'product minimum m/z', args['fminmz'])
    _update(settings, 'search', 'peptide maximum length', args['lmax'])
    _update(settings, 'search', 'peptide minimum length', args['lmin'])
    _update(settings, 'search', 'peptide maximum mass', args['massmax'])
    _update(settings, 'search', 'peptide minimum mass', args['massmin'])
    _update(settings, 'search', 'enzyme', args['e'])
    _update(settings, 'search', 'number of missed cleavages', args['mc'])
    _update(settings, 'search', 'maximum charge', args['cmax'])
    _update(settings, 'search', 'minimum charge', args['cmin'])
    _update(settings, 'search', 'maximum unknown charge', args['cumax'])
    _update(settings, 'search', 'minimum unknown charge', args['cumin'])
    _update(settings, 'search', 'precursor isotope mass error', args['ime'])
    _update(settings, 'search', 'shifts', args['shifts'])
    _update(settings, 'search', 'snp', args['snp'])
    _update(settings, 'output', 'minimum matched', args['mm'])
    _update(settings, 'input', 'add decoy', ('no' if not args['ad'] else 'yes'))
    _update(settings, 'input', 'decoy prefix', args['prefix'])
    _update(settings, 'input', 'decoy method', args['method'])
    _update(settings, 'input', 'deisotope', args['deis'])
    _update(settings, 'input', 'deisotoping mass tolerance', args['deistol'])
    _update(settings, 'scoring', 'score', 'identipy.scoring.' + args['score'])
    _update(settings, 'scoring', 'minimum peaks', args['minp'])
    _update(settings, 'scoring', 'maximum peaks', args['maxp'])
    _update(settings, 'scoring', 'dynamic range', args['dyn'])
    _update(settings, 'scoring', 'maximum fragment charge', args['mfc'])
    _update(settings, 'performance', 'processes', args['nproc'])
    _update(settings, 'modifications', 'maximum variable mods', args['maxmods'])
    _update(settings, 'modifications', 'protein nterm cleavage', args['ncleave'])
    _update(settings, 'modifications', 'protein cterm cleavage', args['ccleave'])

    inputfile = args['file']

    utils.write_output(inputfile, settings, main.process_file(inputfile, settings))
