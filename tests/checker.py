from identipy import main, utils
from sys import argv

settings = main.settings('test.cfg')

try:
    inputfile = argv[1] if len(argv) > 1 else 'test.mgf'
    utils.write_pepxml(inputfile, settings, main.process_file(inputfile, settings))
except KeyboardInterrupt:
    pass
finally:
    print('The search is finished')
