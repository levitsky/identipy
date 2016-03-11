import os
from .spectrum_centric import *
from .peptide_centric import *
from . import utils


def process_file(fname, settings):
    
    stage1 = settings.get('misc', 'first stage')
    if stage1:
        return double_run(fname, settings, utils.import_(stage1))
    else:
        iterate = settings.get('misc', 'iterate')
        ftype = fname.rsplit('.', 1)[-1].lower()
        if iterate == 'spectra':
            spectra = utils.iterate_spectra(fname)
            return process_spectra(spectra, settings)
        elif iterate == 'peptides':
            utils.seen_target.clear()
            utils.seen_decoy.clear()
            return process_peptides(fname, settings)
        else:
            raise ValueError('iterate must be "spectra" or "peptides"')

def double_run(fname, settings, stage1):
    print '[double run] stage 1 starting ...'
    new_settings = stage1(fname, settings)
    print '[double run] stage 2 starting ...'
    return process_file(fname, new_settings)


def settings(fname=None, default_name=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'default.cfg')):
    """Read a configuration file and return a :py:class:`RawConfigParser` object.
    """

    raw_config = utils.CustomRawConfigParser(dict_type=dict, allow_no_value=True)
    if default_name:
        raw_config.read(default_name)
    if fname:
        raw_config.read(fname)
    return raw_config

