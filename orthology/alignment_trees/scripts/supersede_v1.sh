#!/bin/bash
cd /home/pdr/work/phase1b
mv set set_A1v1_superseded && mv aln aln_A1v1_superseded && mv trees trees_A1v1_superseded
mkdir -p set aln trees
cp -r set_A1v1_superseded/sampling set/sampling
ls
