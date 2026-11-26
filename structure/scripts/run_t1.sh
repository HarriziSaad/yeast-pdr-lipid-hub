#!/bin/bash
export MAMBA_ROOT_PREFIX=/home/pdr/micromamba
P=/mnt/c/Users/MSI/Documents/Advisor/project
/home/pdr/bin/micromamba run -n pdr python -I $P/phase3_structure/scripts/t1_download.py $P/phase3_structure/downloads $P/phase0_verify/residues/downloads 2>&1 | tail -20
