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

def main():
    parser = argparse.ArgumentParser(
        description='Search proteins using LC-MS/MS spectra',
        epilog='''

    Example usage
    -------------
    $ identipy input.mgf -db human.fasta
    -------------
    ''',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('file', help='input .mzML or .mgf file with MS/MS spectra')
    parser.add_argument('-db', help='path to protein fasta file', required=True)
    parser.add_argument('-cfg', help='path to file with parameters')
    parser.add_argument('-punit', help='precursor mass tolerance unit. Can be ppm or Da', default='ppm', type=str)
    parser.add_argument('-ptol', help='precursor mass tolerance', default=10.0, type=float)
    parser.add_argument('-lptol', help='*left precursor mass tolerance', default=False, type=float)
    parser.add_argument('-rptol', help='*right precursor mass tolerance', default=False, type=float)
    parser.add_argument('-funit', help='fragment mass tolerance unit. Can be ppm or Da', default='Da')
    parser.add_argument('-ftol', help='fragment mass tolerance', default=0.3, type=float)
    parser.add_argument('-fminmz', help='fragment min m/z', default=100.0, type=float)
    parser.add_argument('-lmin', help='min length of peptides', default=5, type=int)
    parser.add_argument('-lmax', help='max length of peptides', default=55, type=int)
    parser.add_argument('-massmin', help='min mass of peptides', default=300, type=int)
    parser.add_argument('-massmax', help='max mass of peptides', default=4000, type=int)
    parser.add_argument('-e', help='cleavage rule in quotes!. X!Tandem style for cleavage rules', default='[RK]|{P}')
    parser.add_argument('-mc', help='number of missed cleavages', default=0, type=int)
    parser.add_argument('-cmin', help='min precursor charge', default=1, type=int)
    parser.add_argument('-cmax', help='max precursor charge', default=9, type=int)
    parser.add_argument('-cumin', help='min unknown precursor charge', default=0, type=int)
    parser.add_argument('-cumax', help='max unknown precursor charge', default=0, type=int)
    parser.add_argument('-ime', help='precursor isotope mass error. The parent ion\
     mass tolerance is expanded by opening up multiple tolerance windows centered\
      on the given number of 13C isotope peaks for a peptide.', default=0, type=int)
    parser.add_argument('-shifts', help='shifts. example: 0,16.000,23.000,12', default='0')
    parser.add_argument('-snp', help='1 means make SNP changes for ALL peptides', default=0, type=int)
    parser.add_argument('-mm', help='number of minimum matched ions', default=4, type=int)
    parser.add_argument('-ad', help='add decoy', default=0, type=int)
    parser.add_argument('-prefix', help='decoy prefix', default='DECOY_')
    parser.add_argument('-method', help='reverse or random', default='reverse')
    parser.add_argument('-deis', help='use MS/MS deisotoping. yes or no', default='no')
    parser.add_argument('-deistol', help='deisotope mass accuracy', default=0.3, type=float)
    parser.add_argument('-score', help='used score. Can be RNHS, hyperscore or morpheusscore', default='RNHS')
    parser.add_argument('-minp', help='minumum peaks in MS/MS spectra', default=4, type=int)
    parser.add_argument('-maxp', help='maximum peaks in MS/MS spectra', default=100, type=int)
    parser.add_argument('-dyn', help='dynamic range', default=100.0, type=float)
    parser.add_argument('-mfc', help='maximum fragment charge', default=1, type=int)
    parser.add_argument('-nproc', help='number of processes. 0 means auto', default=0, type=int)
    parser.add_argument('-maxmods', help='maximum variable mods per sequence', default=2, type=int)
    parser.add_argument('-ncleave', help='protein nterm cleavage', default=1.007825, type=float)
    parser.add_argument('-ccleave', help='protein cterm cleavage', default=17.002735, type=float)
    parser.add_argument('-fmods', help='fixed modifications. in mass1@aminoacid1,mass2@aminoacid2 format', default='57.021464@C')
    parser.add_argument('-vmods', help='variable modifications. in mass1@aminoacid1,mass2@aminoacid2 format', default='')

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

    settings.set('input', 'database', args['db'])
    settings.set('search', 'precursor accuracy unit', args['punit'])
    settings.set('search', 'precursor accuracy left', (args['ptol'] if not args['lptol'] else args['lptol']))
    settings.set('search', 'precursor accuracy right', (args['ptol'] if not args['rptol'] else args['rptol']))
    settings.set('search', 'product accuracy unit', args['funit'])
    settings.set('search', 'product accuracy', args['ftol'])
    settings.set('search', 'product minimum m/z', args['fminmz'])
    settings.set('search', 'peptide maximum length', args['lmax'])
    settings.set('search', 'peptide minimum length', args['lmin'])
    settings.set('search', 'peptide maximum mass', args['massmax'])
    settings.set('search', 'peptide minimum mass', args['massmin'])
    settings.set('search', 'enzyme', args['e'])
    settings.set('search', 'number of missed cleavages', args['mc'])
    settings.set('search', 'maximum charge', args['cmax'])
    settings.set('search', 'minimum charge', args['cmin'])
    settings.set('search', 'maximum unknown charge', args['cumax'])
    settings.set('search', 'minimum unknown charge', args['cumin'])
    settings.set('search', 'precursor isotope mass error', args['ime'])
    settings.set('search', 'shifts', args['shifts'])
    settings.set('search', 'snp', args['snp'])
    settings.set('output', 'minimum matched', args['mm'])
    settings.set('input', 'add decoy', ('no' if not args['ad'] else 'yes'))
    settings.set('input', 'decoy prefix', args['prefix'])
    settings.set('input', 'decoy method', args['method'])
    settings.set('input', 'deisotope', args['deis'])
    settings.set('input', 'deisotoping mass tolerance', args['deistol'])
    settings.set('scoring', 'score', 'identipy.scoring.' + args['score'])
    settings.set('scoring', 'minimum peaks', args['minp'])
    settings.set('scoring', 'maximum peaks', args['maxp'])
    settings.set('scoring', 'dynamic range', args['dyn'])
    settings.set('scoring', 'maximum fragment charge', args['mfc'])
    settings.set('performance', 'processes', args['nproc'])
    settings.set('modifications', 'maximum variable mods', args['maxmods'])
    settings.set('modifications', 'protein nterm cleavage', args['ncleave'])
    settings.set('modifications', 'protein cterm cleavage', args['ccleave'])

    inputfile = args['file']
# cfg = argv[1]
    utils.write_pepxml(inputfile, settings, main.process_file(inputfile, settings))
