#!/bin/sh

cd test_data
identipy -cfg identipy.cfg -o . *.mgf
pyteomics pepxml info *.pep.xml

echo "Reference values (1% FDR):"
echo 26718
echo 17670
