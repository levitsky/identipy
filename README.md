**IdentiPy** is a search engine for bottom-up proteomics written in Python.

# How to install #

```
#!shell

$ hg clone https://bitbucket.org/levitsky/identipy
$ cd identipy
$ pip install .

```

# Requirements #

 - Python 2.7
 - scipy
 - pyteomics
 - lxml
 
Not strictly required, but highly recommended:

 - cython
 - pyteomics.cythonize

# How to use #

For help on command-line usage, run:

```
#!shell

$ identipy --help
```

Search parameters can be specified using command-line options or in a configuration file.
Allowed parameters and their default values are listed in the
[default configuration file](https://bitbucket.org/levitsky/identipy/src/tip/identipy/default.cfg).


# Related projects #


 - Pyteomics: https://bitbucket.org/levitsky/pyteomics

 - pyteomics.cythonize: https://github.com/mobiusklein/pyteomics.cythonize

 - MP score: https://bitbucket.org/markmipt/mp-score

 - IdentiPy Server: https://bitbucket.org/levitsky/identipy_server