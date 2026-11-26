#!/bin/bash
set -e
O=/mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1b/outputs/set
mkdir -p $O/sampling
cd /home/pdr/work/phase1b/set
cp sampling/sampling_plan.tsv sampling/sampling_summary.json sampling/dias_match.tsv sampling/ygob_match.tsv $O/sampling/
cp analysis_set.faa analysis_set_manifest.tsv analysis_set_summary.json counts_order_x_family.tsv counts_besthit_label.tsv counts_besthit_species.tsv $O/
cd $O
sha256sum analysis_set.faa analysis_set_manifest.tsv sampling/sampling_plan.tsv analysis_set_summary.json
grep -c '>' analysis_set.faa
