#!/bin/sh
rm -rf 0
rm -rf 0.00
rm -rf *[1-9]*
rm -rf processor*[0-9]*
rm -rf probes*[0-9]*
rm -f  log.txt
rm -rf postProcessing
rm -rf dynamicCode
rm -rf constant/polyMesh
rm -f  mass_balance_check.csv
find . -name "*~" -type f -delete
