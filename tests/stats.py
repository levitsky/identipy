from identipy import main, scoring
from sys import argv
import numpy as np

inputfile = argv[1] if len(argv) > 1 else 'test.mgf'
settings = main.settings(inputfile.replace('.mgf', '.cfg'))

for res in main.process_file(inputfile, settings):
    scores = 8 * np.log10(np.array([c[0] for c in res['candidates']]))
    if scores.size: print scores[0]
    coeffs = scoring.survival_hist(scores)[1]
    print coeffs


