import argparse
import string
import logging.config
import os
import subprocess
import copy
import shlex

LOGGING = {
    'version': 1,
    'disable_existing_loggers': True,
    'formatters': {
        'verbose': {
            'format': '%(levelname)7s %(asctime)s %(module)s %(process)d %(thread)d %(message)s',
        },
        'simple': {
            'format': '%(levelname)7s: %(asctime)s %(message)s',
            'datefmt': '[%H:%M:%S]',
        },
    },
    'handlers': {
        'console': {
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'simple',
        },
    },
    'loggers': {
        'identipy': {
            'handlers': ['console'],
            'level': 'INFO',
        }
    }
}

logging.config.dictConfig(LOGGING)
import logging
logger = logging.getLogger(__name__)
from . import main, utils


def get_label(modmass, labels):
    abt = string.ascii_lowercase
    abt_l = len(abt) - 1
    if modmass in labels:
        return labels[modmass], labels, 0
    else:
        labels[modmass] = abt[labels['i']] + abt[labels['j']] + abt[labels['k']]
        labels['k'] += 1
        if labels['k'] > abt_l:
            labels['k'] = 0
            labels['j'] += 1
        if labels['j'] > abt_l:
            labels['j'] = 0
            labels['i'] += 1
        return labels[modmass], labels, 1


def process_mods(settings, spec, name, labels):
    mods_array = []
    if spec:
        for mod in spec.split(','):
            modmass, modaa = mod.split('@')
            lbl, labels, flag = get_label(modmass, labels)
            if modaa == '[':
                ntermlabel, modaa, ctermlabel = '-', '', ''
            elif modaa == ']':
                ntermlabel, modaa, ctermlabel = '', '', '-'
            else:
                ntermlabel, ctermlabel = '', ''
            mods_array.append(ctermlabel + lbl + modaa + ntermlabel)
            if flag:
                settings.set('modifications', lbl, modmass)
    if mods_array or spec is not None:
        settings.set('modifications', name, ','.join(mods_array))


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

    parser.add_argument('file',     help='input mzML or MGF file with MS/MS spectra', nargs='+')
    parser.add_argument('-db',      help='path to protein FASTA file', metavar='FASTA')
    parser.add_argument('-cfg',     help='path to file with parameters', metavar='CONFIG_FILE')
    parser.add_argument('-out', '-o', help='output path', metavar='PATH')
    parser.add_argument('-of',      help='output format', metavar='FORMAT')
    parser.add_argument('-sep',     help='output column separator (for table format)')
    parser.add_argument('-at',      help='Use auto-tuning of search parameters', action='store_true')
    parser.add_argument('-nopwide', help='Do not increase initial precursor mass accuracy for auto-tuning', action='store_true')
    parser.add_argument('-punit',   help='precursor mass tolerance unit', metavar='UNIT', choices=['ppm', 'Da'])
    parser.add_argument('-ptol',    help='precursor mass tolerance', type=float, metavar='VALUE')
    parser.add_argument('-lptol',   help='*left precursor mass tolerance', type=float, metavar='VALUE')
    parser.add_argument('-rptol',   help='*right precursor mass tolerance', type=float, metavar='VALUE')
    parser.add_argument('-funit',   help='fragment mass tolerance unit', metavar='UNIT', choices=['ppm', 'Da'])
    parser.add_argument('-ftol',    help='fragment mass tolerance', type=float, metavar='VALUE')
    parser.add_argument('-fminmz',  help='fragment min m/z', type=float, metavar='VALUE')
    parser.add_argument('-lmin',    help='min length of peptides', type=int, metavar='N')
    parser.add_argument('-lmax',    help='max length of peptides', type=int, metavar='N')
    parser.add_argument('-massmin', help='min mass of peptides', type=float, metavar='VALUE')
    parser.add_argument('-massmax', help='max mass of peptides', type=float, metavar='VALUE')
    parser.add_argument('-e',       help='cleavage rule in quotes!. X!Tandem style for cleavage rules', metavar='RULE')
    parser.add_argument('-mc',      help='number of missed cleavages', type=int, metavar='N')
    parser.add_argument('-semi',    help='include semitryptic peptides', action='store_true')
    parser.add_argument('-noclip',  help='Disable clipping of N-terminal methionine', action='store_false', dest='clip_M')
    parser.add_argument('-cmin',    help='min precursor charge', type=int, metavar='N')
    parser.add_argument('-cmax',    help='max precursor charge', type=int, metavar='N')
    parser.add_argument('-cumin',   help='min unknown precursor charge', type=int, metavar='N')
    parser.add_argument('-cumax',   help='max unknown precursor charge', type=int, metavar='N')
    parser.add_argument('-ime',     help='precursor isotope mass error. The parent ion\
     mass tolerance is expanded by opening up multiple tolerance windows centered\
      on the given number of 13C isotope peaks for a peptide.', type=int, metavar='N')
    parser.add_argument('-shifts',  help='shifts. example: 0,16.000,23.000,12')
    parser.add_argument('-snp',     help='1 means make SNP changes for ALL peptides', type=int)
    parser.add_argument('-rapid',   help='leave only 2000 random spectra for processing', action='store_true')
    parser.add_argument('-mm',      help='number of minimum matched ions', type=int, metavar='N')
    parser.add_argument('-ad',      help='add decoy', action='store_true')
    parser.add_argument('-prefix',  help='decoy prefix')
    parser.add_argument('-infix',   help='decoy infix')
    parser.add_argument('-method',  help='decoy method; reverse or shuffle', choices=['reverse', 'shuffle'])
    parser.add_argument('-nodeis',  help='do not use MS/MS deisotoping', action='store_true')
    parser.add_argument('-deistol', help='deisotope mass accuracy', type=float)
    parser.add_argument('-score',   help='used scoring function', choices=['RNHS2', 'RNHS', 'hyperscore', 'morpheusscore'])
    parser.add_argument('-minp',    help='minumum peaks in MS/MS spectra', type=int, metavar='N')
    parser.add_argument('-maxp',    help='maximum peaks in MS/MS spectra', type=int, metavar='N')
    parser.add_argument('-dyn',     help='dynamic range', type=float)
    parser.add_argument('-mfc',     help='maximum fragment charge', type=int, metavar='N')
    parser.add_argument('-nproc',   help='number of processes. 0 means auto', type=int, metavar='N')
    parser.add_argument('-maxmods', help='maximum variable mods per sequence', type=int, metavar='N')
    parser.add_argument('-ncleave', help='protein nterm cleavage', type=float)
    parser.add_argument('-ccleave', help='protein cterm cleavage', type=float)
    parser.add_argument('-fmods',   help='fixed modifications. Format: mass1@aminoacid1,mass2@aminoacid2')
    parser.add_argument('-vmods',   help='variable modifications. Format: mass1@aminoacid1,mass2@aminoacid2')
    parser.add_argument('-pmods',   help='variable protein terminal modifications')
    parser.add_argument('-tags',    help='Add quantitation tags to the pepXML output. Can be tmt10plex, tmt6plex, tmt11plex or custom format label1:mass1,label2:mass2...')
    parser.add_argument('-debug',   help='Print debugging messages', action='store_true')
    parser.add_argument('-dino',    help='path to Dinosaur JAR file or Biosaur executable. Used for chimeric spectrum processing and MS1 Intensity calculation', default=False)
    parser.add_argument('-dinoargs', help='extra arguments to Dinosaur or Biosaur.', default='')
    parser.add_argument('-sd', '-skipdino', action='store_true', help='Skip feature detection if a feature file is found.')
    parser.add_argument('-demixing',help='Use demixing', action='store_true')
    parser.add_argument('-pif',     help='Calculate PIF', action='store_true')

    args = vars(parser.parse_args())
    if args['debug']:
        logging.getLogger('identipy').setLevel(logging.DEBUG)

    if args['cfg']:
        settings = main.settings(args['cfg'])
    else:
        settings = main.settings()

    labels = {'i': 0, 'j': 0, 'k': 0}
    process_mods(settings, args['fmods'], 'fixed', labels)
    process_mods(settings, args['vmods'], 'variable', labels)
    process_mods(settings, args['pmods'], 'protein variable', labels)

    _update(settings, 'input',  'database', args['db'])
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
    _update(settings, 'search', 'semitryptic', args['semi'])
    _update(settings, 'search', 'clip N-terminal methionine', str(args['clip_M']))
    _update(settings, 'search', 'maximum charge', args['cmax'])
    _update(settings, 'search', 'minimum charge', args['cmin'])
    _update(settings, 'search', 'maximum unknown charge', args['cumax'])
    _update(settings, 'search', 'minimum unknown charge', args['cumin'])
    _update(settings, 'search', 'precursor isotope mass error', args['ime'])
    _update(settings, 'search', 'shifts', args['shifts'])
    _update(settings, 'search', 'snp', args['snp'])
    _update(settings, 'output', 'minimum matched', args['mm'])
    if args['ad']:
        _update(settings, 'input', 'add decoy', 'yes')
    if args['rapid']:
        _update(settings, 'search', 'rapid_check', 1)
    _update(settings, 'input',  'decoy prefix', args['prefix'])
    _update(settings, 'input',  'decoy infix', args['infix'])
    _update(settings, 'input',  'decoy method', args['method'])
    if args['nodeis']:
        _update(settings, 'input',  'deisotope', 'no')
    _update(settings, 'input',  'deisotoping mass tolerance', args['deistol'])
    if args['score']:
        _update(settings, 'scoring', 'score', 'identipy.scoring.' + args['score'])
    _update(settings, 'scoring', 'minimum peaks', args['minp'])
    _update(settings, 'scoring', 'maximum peaks', args['maxp'])
    _update(settings, 'scoring', 'dynamic range', args['dyn'])
    _update(settings, 'scoring', 'maximum fragment charge', args['mfc'])
    _update(settings, 'performance', 'processes', args['nproc'])
    _update(settings, 'modifications', 'maximum variable mods', args['maxmods'])
    _update(settings, 'modifications', 'protein nterm cleavage', args['ncleave'])
    _update(settings, 'modifications', 'protein cterm cleavage', args['ccleave'])
    _update(settings, 'output', 'path', args['out'])
    _update(settings, 'output', 'format', args['of'])
    _update(settings, 'output', 'separator', args['sep'])
    _update(settings, 'output', 'tags', args['tags'])
    if args['at']:
        ao_setting = 'identipy.extras.optimization'
        if args['nopwide']:
            _update(settings, 'optimization', 'increase precursor mass tolerance', 'no')
    else:
        ao_setting = None
    _update(settings, 'misc', 'first stage', ao_setting)

    dino_path = args['dino']
    demixing = args['demixing']
    calc_PIF = args['pif']

    for inputfile in args['file']:
        csettings = copy.deepcopy(settings)

        if dino_path or calc_PIF:
            logger.info('Starting mzML analysis...')
            if os.path.splitext(inputfile)[1].lower() != '.mzml':
                if dino_path:
                    logger.error('Only mzML supported for Dinosaur!')
                elif calc_PIF:
                    logger.error('mzML required for PIF calculation!')
            else:
                try:
                    if dino_path:
                        path_to_features = os.path.splitext(inputfile)[0] + os.extsep + 'features' + os.extsep + 'tsv'
                        if not args['skipdino'] or not os.path.exists(path_to_features):
                            if dino_path.endswith('.jar'):
                                advpath = '--advParams=' + os.path.join(os.path.dirname(os.path.realpath(__file__)), 'adv.txt')
                                logger.info('Starting Dinosaur...')
                                subprocess.call(['java', '-Djava.awt.headless=true', '-jar', os.path.realpath(dino_path), advpath, '--concurrency=12', inputfile] + args['dinoargs'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            elif 'dinosaur' in dino_path:
                                advpath = '--advParams=' + os.path.join(os.path.dirname(os.path.realpath(__file__)), 'adv.txt')
                                logger.info('Starting Dinosaur...')
                                subprocess.call([os.path.realpath(dino_path), advpath, '--concurrency=12', inputfile] + shlex.split(args['dinoargs']), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            elif 'biosaur2' in dino_path:
                                logger.info('Starting biosaur2...')
                                cmd = [os.path.realpath(dino_path), inputfile, '-o', path_to_features] + shlex.split(args['dinoargs'])
                                logger.debug('Running command: %s', cmd)
                                subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            else:
                                logger.info('Starting Biosaur...')
                                subprocess.call([os.path.realpath(dino_path), inputfile, '-out', path_to_features] + shlex.split(args['dinoargs']), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        if demixing:
                            logger.info('Starting demultiplexing...')
                    else:
                        path_to_features = False
                    path_to_mgf = utils.demix_chimeric(path_to_features, inputfile, demixing, calc_PIF)
                    logger.info('MGF was created.')
                    if demixing:
                        logger.info('Demultiplexing has finished.')
                    utils.write_output(path_to_mgf, csettings, main.process_file(path_to_mgf, csettings))
                    return
                except Exception as e:
                    logger.error(e)

        utils.write_output(inputfile, csettings, main.process_file(inputfile, csettings))


if __name__ == '__main__':
    run()
