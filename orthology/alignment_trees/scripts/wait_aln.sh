#!/bin/bash
log=/home/pdr/work/phase1b/logs/$1.stderr.log
until grep -q '^exit=' $log 2>/dev/null; do sleep 20; done
grep '^exit=' $log
ls -la /home/pdr/work/phase1b/aln
