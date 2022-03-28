#!/usr/bin/env python

'''
setup.py file for identipy
'''
import os
from setuptools import setup, Extension
import subprocess
import sys

def get_version():
    try:
        version = subprocess.check_output(['git', 'describe']).strip().decode('ascii').replace('-', '.')
        if version[0] == 'v':
            version = version[1:]
        head, tail = version.rsplit('.', 1)
        if not tail.isdigit():
            version = head
    except subprocess.CalledProcessError:
        version = open('VERSION').readline().strip()
    return version


def make_extensions():
    is_ci = bool(os.getenv("CI", ""))
    include_diagnostics = False
    try:
        import numpy
    except ImportError:
        print("C Extensions require `numpy`")
        raise
    try:
        from pyteomics import _capi
    except ImportError:
        print("C Extensions require `pyteomics.cythonize`")
        raise
    try:
        from Cython.Build import cythonize
        cython_directives = {
            'embedsignature': True,
            'profile': include_diagnostics,
            'language_level': sys.version_info.major
        }
        macros = []
        if include_diagnostics:
            macros.append(("CYTHON_TRACE_NOGIL", "1"))
        if is_ci and include_diagnostics:
            cython_directives['linetrace'] = True

        extensions = cythonize([
            Extension(name='identipy.cparser', sources=['identipy/cparser.pyx']),
            Extension(name='identipy.cutils', sources=['identipy/cutils.pyx'],
                      include_dirs=[numpy.get_include(), _capi.get_include()])
        ], compiler_directives=cython_directives)
    except ImportError:
        extensions = [
            Extension(name='identipy.cparser', sources=['identipy/cparser.c']),
            Extension(name='identipy.cutils', sources=['identipy/cutils.c'],
                      include_dirs=[numpy.get_include(), _capi.get_include()])

        ]
    return extensions


def do_setup(cext=True):
    setup(
        name             = 'identipy',
        version          = get_version(),
        description      = '''Pyteomics-based search engine''',
        long_description = (''.join(open('README.md').readlines()) + '\n'
                            + ''.join(open('INSTALL').readlines())),
        author           = 'Lev Levitsky & Mark Ivanov',
        author_email     = 'pyteomics@googlegroups.com',
        url              = 'https://github.com/levitsky/identipy',
        packages         = ['identipy', ],
        package_data     = {'identipy': ['default.cfg', ]},
        install_requires = [line.strip() for line in open('requirements.txt')],
        ext_modules      = make_extensions() if cext else None,
        classifiers      = ['Intended Audience :: Science/Research',
                            'Programming Language :: Python :: 2.7',
                            'Topic :: Scientific/Engineering :: Bio-Informatics',
                            'Topic :: Software Development :: Libraries'],
        license          = 'License :: OSI Approved :: Apache Software License',
        entry_points     = {'console_scripts': ['identipy = identipy.cli:run',
                                                'identipy2pin = identipy.identipy2pin:run']}
        )


try:
    do_setup(True)
except Exception as err:
    print("*" * 60)
    print("Could not compile C Extensions due to %r, attempting pure Python installation." % (err,))
    print("*" * 60)
    do_setup(False)
    print("Could not compile C Extensions due to %r, speedups are not enabled." % (err,))
