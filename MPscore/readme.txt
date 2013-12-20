Script is written in python version 2.7. This script require numpy, scipy, matplotlib, pyteomics and pyteomics.biolccc modules.

The .pep.xml, .fasta and .cfg files are required for basic work of the script. Cfg file contains settings for the algorithm. For some descriptors .mgf and .mzml files are needed. For quantitation through sum of matched fragment ions intensities .t.xml X! Tandem output file is needed.

Algorithm can be run with following command: python MP.py path_to_pepxml path_to_fasta path_to_cfg *path_to_txml *path_to_mzml *path_to_mgf.

Output contains pictures with descriptors distributions, .csv table with identified proteins groups (only best protein from the group in terms of number of identified peptides, at least 2 peptides) and .csv table with identified peptides (only best match for peptide in terms of e-value).
