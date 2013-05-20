from identipy import main as ip
from identipy.scoring import hyperscore
from pyteomics import mgf
import sys
import cProfile
if len(sys.argv) > 1:
    num = int(sys.argv[1])
else:
    num = 3
settings = ip.settings('test.cfg')
cProfile.run(
    '''for _, (s, res) in zip(range(num), ip.process_file('swedcad.mgf', settings)):
        pass
    ''')
