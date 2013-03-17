from identipy import main as ip
from identipy.scoring import hyperscore
from pyteomics import mgf
import sys
import cProfile
if len(sys.argv) > 1:
    num = int(sys.argv[1])
with mgf.read('swedcad.mgf') as r:
    cProfile.run(
        '''for _, (s, r) in zip(range(num), ip.process_file(r, hyperscore, 1)):
            pass
        ''')
