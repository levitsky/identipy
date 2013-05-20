from __future__ import print_function
from pyteomics import mgf
from identipy import main, scoring
settings = main.settings('test.cfg')
true = total = 0
try:
    for result in main.process_file('swedcad.mgf', settings):
        total += 1
        if result['candidates'] and (result['candidates'][0][1] == result['spectrum']['params']['title']):
            true += 1
        else:
            print('{} != {}'.format(result['candidates'][0][1] if result['candidates'] else 'NOTHING', result['spectrum']['params']['title']))
except KeyboardInterrupt:
    pass
finally:
    print('\rCorrect: {} of {}'.format(true, total))

