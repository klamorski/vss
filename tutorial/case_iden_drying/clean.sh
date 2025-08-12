#!/bin/sh
rm -rf *[1-9]*
rm -rf processor*[0-9]*
rm -rf probes*[0-9]*
rm -f  log.txt
rm -rf postProcessing
rm -rf dynamicCode
find . -name "*~" -type f -delete
