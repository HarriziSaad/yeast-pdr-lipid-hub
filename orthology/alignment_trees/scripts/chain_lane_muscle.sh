#!/bin/bash
until grep -q 'done fast+G1' /home/pdr/work/phase1b/logs/model_muscle.txt 2>/dev/null; do sleep 30; done
S=/mnt/c/Users/MSI/Documents/Advisor/project/phase1_orthology/phase1b/scripts
bash $S/full_lane.sh muscle 6 > /home/pdr/work/phase1b/logs/lane_muscle.log 2>&1
