#!/bin/bash
snakemake --configfile m5U_Detector.yaml -s m5U_Detector.py --dag | dot -Tpdf > dag.pdf

#snakemake --dag report.html | dot -Tsvg > dag.svg
#snakemake -s *** --cores 10
#snakemake --cores 10(-j 10)

