import os
from . import utils, peptide_centric
from pyteomics import fasta
import tempfile
import logging
import re
logger = logging.getLogger(__name__)

def process_file(fname, settings):

    fmods = settings.get('modifications', 'fixed')
    if fmods:
        newfixed = []
        for mod in re.split(r'[,;]\s*', fmods):
            if mod.startswith('-'):
                mod_label = mod[1:]
                mass_change = settings.getfloat('modifications', mod_label)
                prev_cterm_mass = settings.getfloat('modifications', 'protein cterm cleavage')
                settings.set('modifications', 'protein cterm cleavage', prev_cterm_mass + mass_change)
            elif mod.endswith('-'):
                mod_label = mod[:-1]
                mass_change = settings.getfloat('modifications', mod_label)
                prev_nterm_mass = settings.getfloat('modifications', 'protein nterm cleavage')
                settings.set('modifications', 'protein nterm cleavage', prev_nterm_mass + mass_change)
            else:
                newfixed.append(mod)
    # settings.set('modifications', 'fixed', ','.join(newfixed))
            # m, aa = parser._split_label(mod)
            # aa_mass[aa] += settings.getfloat('modifications', m)

    add_decoy = settings.getboolean('input', 'add decoy')
    prefix = settings.get('input', 'decoy prefix')
    mode = settings.get('input', 'decoy method')
    db = settings.get('input', 'database')
    if add_decoy and utils.is_db_target_only(db, prefix):
        ft = tempfile.NamedTemporaryFile(mode='w', delete=False)
        fasta.write_decoy_db(db, ft, mode=mode, prefix=prefix)
        ft.flush()
        settings.set('input', 'database', ft.name)
        settings.set('input', 'add decoy', 'no')
        logger.debug('Temporary database: %s (%s)', ft.name, os.path.isfile(ft.name))
    stage1 = settings.get('misc', 'first stage')
    if stage1:
        return double_run(fname, settings, utils.import_(stage1))
    else:
        utils.seen_target.clear()
        utils.seen_decoy.clear()
        return peptide_centric.process_peptides(fname, settings)
        

def double_run(fname, settings, stage1):
    logger.info('[double run] stage 1 starting ...')
    settings.set('misc', 'fast first stage', 1)
    new_settings = stage1(fname, settings)
    logger.info('[double run] stage 2 starting ...')
    new_settings.set('misc', 'fast first stage', 0)
    return process_file(fname, new_settings)


def settings(fname=None, default_name=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'default.cfg')):
    """Read a configuration file and return a :py:class:`RawConfigParser` object.
    """
    raw_config = utils.CustomRawConfigParser(dict_type=dict, allow_no_value=True)
    if default_name:
        logger.info('Reading defaults from %s', default_name)
        if not os.path.isfile(default_name):
            logger.error('FILE NOT FOUND: %s', default_name)
        raw_config.read(default_name)
    if fname:
        logger.info('Reading config from %s', fname)
        if not os.path.isfile(fname):
            logger.error('FILE NOT FOUND: %s', fname)
        raw_config.read(fname)

    acc_unit = raw_config.get('search', 'product accuracy unit')
    if acc_unit == 'ppm':
        acc_ppm = raw_config.getfloat('search', 'product accuracy')
        acc_raw = acc_ppm / 1e6 * 2000
        raw_config.set('search', 'product accuracy', acc_raw)
        raw_config.set('search', 'product accuracy ppm', acc_ppm)

    return raw_config
