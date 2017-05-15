from identipy import main, utils
from sys import argv

inputfile = argv[2]
cfg = argv[1]
settings = main.settings(cfg)
utils.write_output(inputfile, settings, main.process_file(inputfile, settings))
