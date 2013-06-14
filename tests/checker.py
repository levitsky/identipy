from pyteomics import parser, fasta
from identipy import main, utils

settings = main.settings('/home/mark/workspace/IdentiPy_new_stableRT/mark.cfg')

aa_mass = utils.get_aa_mass(settings)

db = settings.get('input', 'database')
enzyme = settings.get('search', 'enzyme')
enzyme = parser.expasy_rules.get(enzyme, enzyme)
mc = settings.getint('search', 'miscleavages')
minlen = settings.getint('search', 'peptide minimum length')
maxlen = settings.getint('search', 'peptide maximum length')

pept_prot = dict()
for desc, prot in fasta.read(db):
    for pep in parser.cleave(prot, enzyme, mc):
        if minlen <= len(pep) <= maxlen and all(aminoacid not in ['B', 'X', 'J', 'Z', 'U', 'O'] for aminoacid in pep):
            dbinfo, note = desc, ('d' if desc.split('|')[0].startswith('DECOY_') else 't')
            try:
                pept_prot[pep].append((dbinfo, note))
            except:
                pept_prot[pep] = [(dbinfo, note)]

try:
    inputfile = '/home/mark/workspace/IdentiPy_new_stableRT/tests/test.mgf'
    utils.write_pepxml(inputfile, settings, main.process_file(inputfile, settings), pept_prot)
except KeyboardInterrupt:
    pass
finally:
    print('The search is finished')
