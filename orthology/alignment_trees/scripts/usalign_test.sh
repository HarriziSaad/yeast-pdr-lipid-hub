#!/bin/bash
mkdir -p /home/pdr/work/phase1b/struct; cd /home/pdr/work/phase1b/struct
P=/mnt/c/Users/MSI/Documents/Advisor/project
cp $P/phase0_verify/residues/downloads/7P04.cif .
cp $P/phase3_structure/downloads/rcsb_cdr1/9VT1.cif .
time USalign 9VT1.cif 7P04.cif -mol prot -outfmt 0 > us_9VT1_7P04.txt 2> us.err
head -40 us_9VT1_7P04.txt | cut -c1-200
tail -c 300 us.err
