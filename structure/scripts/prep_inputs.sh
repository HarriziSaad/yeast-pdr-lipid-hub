#!/bin/bash
set -e
P=/mnt/c/Users/MSI/Documents/Advisor/project
W=/home/pdr/work/phase3
mkdir -p $W/in $W/scripts $W/out
rsync -a $P/phase3_structure/downloads/ $W/in/
rsync -a $P/phase3_structure/scripts/ $W/scripts/
cp $P/phase0_verify/residues/scripts/harris_th_helices.tsv $W/in/
cp $P/preregistration/site_classes_FROZEN.tsv $W/in/
cp $P/preregistration/outputs/structure_covariates_7P04.tsv $W/in/
ls $W/in | head
cp $P/preregistration/downloads/opm_7p04/opm_7p04.pdb $W/in/
