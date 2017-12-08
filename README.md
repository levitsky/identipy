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

## GUI way ##

You can separately install a web-based GUI for IdentiPy, [IdentiPy Server](https://bitbucket.org/levitsky/identipy_server).
Please refer to the linked page for system requirements and installation instructions.

## CLI way ##

A typical command to process a file would look like this:

```
#!shell

$ identipy -cfg my.cfg spectra.mgf
```
Here, `my.cfg` is a settings file specifying the search parameters. Allowed parameters and their default values are listed in the
[default configuration file](https://bitbucket.org/levitsky/identipy/src/tip/identipy/default.cfg).
Settings not specified in `my.cfg` will be taken from the default file.

Search settings can also be overriden using command-line options.

For help on command-line usage, run:

```
#!shell

$ identipy --help
```

You will see a message like this:

```

$ identipy --help
usage: identipy [-h] [-db DB] [-cfg CFG] [-punit PUNIT] [-ptol PTOL]
                [-lptol LPTOL] [-rptol RPTOL] [-funit FUNIT] [-ftol FTOL]
                [-fminmz FMINMZ] [-lmin LMIN] [-lmax LMAX] [-massmin MASSMIN]
                [-massmax MASSMAX] [-e E] [-mc MC] [-cmin CMIN] [-cmax CMAX]
                [-cumin CUMIN] [-cumax CUMAX] [-ime IME] [-shifts SHIFTS]
                [-snp SNP] [-mm MM] [-ad AD] [-prefix PREFIX] [-method METHOD]
                [-deis DEIS] [-deistol DEISTOL] [-score SCORE] [-minp MINP]
                [-maxp MAXP] [-dyn DYN] [-mfc MFC] [-nproc NPROC]
                [-maxmods MAXMODS] [-ncleave NCLEAVE] [-ccleave CCLEAVE]
                [-fmods FMODS] [-vmods VMODS]
                file

Search proteins using LC-MS/MS spectra

positional arguments:
  file              input .mzML or .mgf file with MS/MS spectra

optional arguments:
  -h, --help        show this help message and exit
  -db DB            path to protein fasta file
  -cfg CFG          path to file with parameters
  -punit PUNIT      precursor mass tolerance unit. Can be ppm or Da
  -ptol PTOL        precursor mass tolerance
  -lptol LPTOL      *left precursor mass tolerance
  -rptol RPTOL      *right precursor mass tolerance
  -funit FUNIT      fragment mass tolerance unit. Can be ppm or Da
  -ftol FTOL        fragment mass tolerance
  -fminmz FMINMZ    fragment min m/z
  -lmin LMIN        min length of peptides
  -lmax LMAX        max length of peptides
  -massmin MASSMIN  min mass of peptides
  -massmax MASSMAX  max mass of peptides
  -e E              cleavage rule in quotes!. X!Tandem style for cleavage
                    rules
  -mc MC            number of missed cleavages
  -cmin CMIN        min precursor charge
  -cmax CMAX        max precursor charge
  -cumin CUMIN      min unknown precursor charge
  -cumax CUMAX      max unknown precursor charge
  -ime IME          precursor isotope mass error. The parent ion mass
                    tolerance is expanded by opening up multiple tolerance
                    windows centered on the given number of 13C isotope peaks
                    for a peptide.
  -shifts SHIFTS    shifts. example: 0,16.000,23.000,12
  -snp SNP          1 means make SNP changes for ALL peptides
  -mm MM            number of minimum matched ions
  -ad AD            add decoy
  -prefix PREFIX    decoy prefix
  -method METHOD    reverse or random
  -deis DEIS        use MS/MS deisotoping. yes or no
  -deistol DEISTOL  deisotope mass accuracy
  -score SCORE      used score. Can be RNHS, hyperscore or morpheusscore
  -minp MINP        minumum peaks in MS/MS spectra
  -maxp MAXP        maximum peaks in MS/MS spectra
  -dyn DYN          dynamic range
  -mfc MFC          maximum fragment charge
  -nproc NPROC      number of processes. 0 means auto
  -maxmods MAXMODS  maximum variable mods per sequence
  -ncleave NCLEAVE  protein nterm cleavage
  -ccleave CCLEAVE  protein cterm cleavage
  -fmods FMODS      fixed modifications. in mass1@aminoacid1,mass2@aminoacid2
                    format
  -vmods VMODS      variable modifications. in
                    mass1@aminoacid1,mass2@aminoacid2 format

    Example usage
    -------------
    $ identipy input.mgf -db human.fasta
    -------------

```

Search parameters can be specified using command-line options or in a configuration file.


# Related projects #


 - Pyteomics: https://bitbucket.org/levitsky/pyteomics

 - pyteomics.cythonize: https://github.com/mobiusklein/pyteomics.cythonize

 - MP score: https://bitbucket.org/markmipt/mp-score

 - IdentiPy Server: https://bitbucket.org/levitsky/identipy_server