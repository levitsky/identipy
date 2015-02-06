#!/usr/bin/env python

'''
setup.py file for identipy
'''

from distutils.core import setup

version = open('VERSION').readline().strip()

setup(
    name = 'identipy',
    version = version,
    description      = '''Pyteomics-based search engine''',
    long_description = (''.join(open('README').readlines()) + '\n'
                        + ''.join(open('INSTALL').readlines())),
    author           = 'Lev Levitsky & Mark Ivanov',
    author_email     = 'pyteomics@googlegroups.com',
    url              = 'http://hg.theorchromo.ru/identipy',
    packages         = ['identipy', ],
    package_data     = {'identipy': ['default.cfg']},
    requires         = [line.strip() for line in open('requirements.txt')],
    classifiers      = ['Intended Audience :: Science/Research',
                        'Programming Language :: Python :: 2.7',
                        'Topic :: Scientific/Engineering :: Bio-Informatics',
                        'Topic :: Software Development :: Libraries'],
    license          = 'License :: OSI Approved :: Apache Software License',
    )
