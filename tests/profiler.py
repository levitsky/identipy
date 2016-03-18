from identipy import main as ip
import cProfile

settings = ip.settings('UPS3.cfg')
if settings.getint('performance', 'processes') != 1:
    settings.set('performance', 'processes', 1)
    print "Setting the number of processes to 1."
cProfile.run(
    '''r = list(ip.process_file('test.mgf', settings))
    ''', 'profiler4.out')
