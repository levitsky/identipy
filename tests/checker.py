from identipy import main, utils
from sys import argv

settings = main.settings('/home/mark/work/IdentiPy/mark.cfg')

try:
    inputfile = argv[1]
    utils.write_pepxml(inputfile, settings, main.process_file(inputfile, settings))
except KeyboardInterrupt:
    pass
finally:
    print('The search is finished')
