from identipy import main, scoring
from sys import argv
import numpy as np
import pylab

inputfile = argv[1] if len(argv) > 1 else 'test.mgf'
settings = main.settings(inputfile.replace('.mgf', '.cfg'))

coeffs = []
for res in main.process_file(inputfile, settings):
    scores = 4 * np.log10(np.array([c[0] for c in res['candidates']]))
#   if scores.size: print 'Top score (log):', scores[0]
    coeffs.append(scoring.survival_hist(scores)[1])

coeffs = np.array(coeffs)
pylab.figure()
pylab.title('Distribution of "a"')
pylab.hist(coeffs[:,0])
pylab.figure()
pylab.title('Distribution of "b"')
pylab.hist(coeffs[:,1])
pylab.show()
