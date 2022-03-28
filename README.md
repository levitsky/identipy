**IdentiPy** is a search engine for bottom-up proteomics written in Python.

# Citation #

IdentiPy is described in this JPR paper: http://dx.doi.org/10.1021/acs.jproteome.7b00640

Please cite it when using IdentiPy or its parts.

# License #

IdentiPy is published under the Apache 2.0 license.

# How to install #

```
pip install git+https://github.com/levitsky/identipy.git
```
or:

```
$ git clone https://github.com/levitsky/identipy
$ cd identipy
$ pip install .

```

# Requirements #

See [requirements.txt](requirements.txt). Key dependencies are:

 - Python
 - scipy
 - pyteomics
 - lxml
 - cython
 - pyteomics.cythonize

# How to use #

## GUI way ##

You can separately install a web-based GUI for IdentiPy, [IdentiPy Server](https://github.com/levitsky/identipy_server).
Please refer to the linked page for system requirements and installation instructions.

## CLI way ##

A typical command to process a file would look like this:

```
$ identipy -cfg my.cfg spectra.mgf
```

Here, `my.cfg` is a settings file specifying the search parameters. Allowed parameters and their default values are listed in the
[default configuration file](identipy/default.cfg).
Settings not specified in `my.cfg` will be taken from the default file.

Search settings can also be overriden using command-line options.

For help on command-line usage, run:

```
$ identipy --help
```

You will see a message like this:

```
$ identipy --help
usage: identipy [-h] [-db FASTA] [-cfg CONFIG_FILE] [-out PATH] [-of FORMAT]
                [-sep SEP] [-at] [-nopwide] [-punit UNIT] [-ptol VALUE]
                [-lptol VALUE] [-rptol VALUE] [-funit UNIT] [-ftol VALUE]
                [-fminmz VALUE] [-lmin N] [-lmax N] [-massmin VALUE]
                [-massmax VALUE] [-e RULE] [-mc N] [-semi] [-noclip] [-cmin N]
                [-cmax N] [-cumin N] [-cumax N] [-ime N] [-shifts SHIFTS]
                [-snp SNP] [-rapid] [-mm N] [-ad] [-prefix PREFIX]
                [-infix INFIX] [-method {reverse,shuffle}] [-deis]
                [-deistol DEISTOL]
                [-score {RNHS2,RNHS,hyperscore,morpheusscore}] [-minp N]
                [-maxp N] [-dyn DYN] [-mfc N] [-nproc N] [-maxmods N]
                [-ncleave NCLEAVE] [-ccleave CCLEAVE] [-fmods FMODS]
                [-vmods VMODS] [-pmods PMODS] [-tags TAGS] [-debug]
                [-dino DINO] [-dinoargs [DINOARGS ...]] [-demixing] [-pif]
                file

Search proteins using LC-MS/MS spectra

positional arguments:
  file                  input .mzML or .mgf file with MS/MS spectra

options:
  -h, --help            show this help message and exit
  -db FASTA             path to protein FASTA file
  -cfg CONFIG_FILE      path to file with parameters
  -out PATH, -o PATH    output path
  -of FORMAT            output format
  -sep SEP              output column separator (for table format)
  -at                   Use auto-tuning of search parameters
  -nopwide              Do not increase initial precursor mass accuracy for
                        auto-tuning
  -punit UNIT           precursor mass tolerance unit
  -ptol VALUE           precursor mass tolerance
  -lptol VALUE          *left precursor mass tolerance
  -rptol VALUE          *right precursor mass tolerance
  -funit UNIT           fragment mass tolerance unit
  -ftol VALUE           fragment mass tolerance
  -fminmz VALUE         fragment min m/z
  -lmin N               min length of peptides
  -lmax N               max length of peptides
  -massmin VALUE        min mass of peptides
  -massmax VALUE        max mass of peptides
  -e RULE               cleavage rule in quotes!. X!Tandem style for cleavage
                        rules
  -mc N                 number of missed cleavages
  -semi                 include semitryptic peptides
  -noclip               Disable clipping of N-terminal methionine
  -cmin N               min precursor charge
  -cmax N               max precursor charge
  -cumin N              min unknown precursor charge
  -cumax N              max unknown precursor charge
  -ime N                precursor isotope mass error. The parent ion mass
                        tolerance is expanded by opening up multiple tolerance
                        windows centered on the given number of 13C isotope
                        peaks for a peptide.
  -shifts SHIFTS        shifts. example: 0,16.000,23.000,12
  -snp SNP              1 means make SNP changes for ALL peptides
  -rapid                leave only 2000 random spectra for processing
  -mm N                 number of minimum matched ions
  -ad                   add decoy
  -prefix PREFIX        decoy prefix
  -infix INFIX          decoy infix
  -method {reverse,shuffle}
                        decoy method; reverse or shuffle
  -deis                 use MS/MS deisotoping
  -deistol DEISTOL      deisotope mass accuracy
  -score {RNHS2,RNHS,hyperscore,morpheusscore}
                        used scoring function
  -minp N               minumum peaks in MS/MS spectra
  -maxp N               maximum peaks in MS/MS spectra
  -dyn DYN              dynamic range
  -mfc N                maximum fragment charge
  -nproc N              number of processes. 0 means auto
  -maxmods N            maximum variable mods per sequence
  -ncleave NCLEAVE      protein nterm cleavage
  -ccleave CCLEAVE      protein cterm cleavage
  -fmods FMODS          fixed modifications. Format:
                        mass1@aminoacid1,mass2@aminoacid2
  -vmods VMODS          variable modifications. Format:
                        mass1@aminoacid1,mass2@aminoacid2
  -pmods PMODS          variable protein terminal modifications
  -tags TAGS            Add quantitation tags to the pepXML output. Can be
                        tmt10plex, tmt6plex, tmt11plex or custom format
                        label1:mass1,label2:mass2...
  -debug                Print debugging messages
  -dino DINO            path to Dinosaur JAR file or Biosaur executable. Used
                        for chimeric spectrum processing and MS1 Intensity
                        calculation
  -dinoargs [DINOARGS ...]
                        extra arguments to Dinosaur or Biosaur.
  -demixing             Use demixing
  -pif                  Calculate PIF

    Example usage
    -------------
    $ identipy input.mgf -db human.fasta
    -------------

```


# Related projects #

 - Pyteomics: https://github.com/levitsky/pyteomics

 - pyteomics.cythonize: https://github.com/mobiusklein/pyteomics.cythonize

 - Scavager: https://github.com/markmipt/scavager

 - IdentiPy Server: https://github.com/levitsky/identipy_server
