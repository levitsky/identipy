#!/bin/sh

cd test_data
find . -name '*.mgf' -print -exec identipy -cfg identipy.cfg -o . {} \;
pyteomics pepxml info *.pep.xml

echo "Reference values (1% FDR):"
echo 26718
echo 17670
