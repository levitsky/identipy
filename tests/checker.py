from identipy import main, utils
from sys import argv

inputfile = argv[1] if len(argv) > 1 else 'test.mgf'
if inputfile == 'swedcad.mgf':
    settings = main.settings('swedcad.cfg')
    tot, wro, leu = 0, 0, 0
    try:
        for res in main.process_file(inputfile, settings):
#           print len(res['candidates'])
            tot += 1
            pep = res['spectrum']['params']['title']
            match = res['candidates'][0][1] if res['candidates'].size else '<NOTHING>'
            if pep != match:
                wro += 1
                if pep.replace('I', 'L') == match.replace('I', 'L'):
                    leu += 1
                else:
                    print pep, '<=>', match
    except KeyboardInterrupt:
        pass
    finally:
        print 'Search finished.'
        if tot:
            print 'FDR: {:.0f}% ({:.0f}% of them Leu/Ile swaps).'.format(
                    100.0 * wro / tot, 100.0 * leu / tot)
        else:
            print 'No results.'

else:
    settings = main.settings('test.cfg')
    utils.write_pepxml(inputfile, settings, main.process_file(inputfile, settings))
    print 'The search is finished.'
