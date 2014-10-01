from identipy import main as ip
from identipy.scoring import hyperscore
from pyteomics import mgf
import sys
import cProfile
if len(sys.argv) > 1:
    num = int(sys.argv[1])
else:
    num = 3
settings = ip.settings('UPS2.cfg')
if settings.getint('performance', 'processes') != 1:
    settings.set('performance', 'processes', 1)
    print "Setting the number of processes to 1."
cProfile.run(
    '''for _, res in zip(range(num), ip.process_file('test.mgf', settings)):
        pass
    ''', 'profiler.out')
