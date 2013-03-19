from __future__ import print_function
from pyteomics import mgf
from identipy import main, scoring
settings = main.settings('test.cfg')
with mgf.read('swedcad.mgf') as r:
    true = total = 0
    try:
        for spectrum, result in main.process_file(r, settings):
            total += 1
            if result and (result[0][1] == spectrum['params']['title']):
                true += 1
            else:
                print('{} != {}'.format(result[0][1], spectrum['params']['title']))
    except KeyboardInterrupt:
        pass
    finally:
        print('\rCorrect: {} of {}'.format(true, total))

