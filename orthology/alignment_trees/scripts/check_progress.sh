#!/bin/bash
cd /home/pdr/work/phase1b
echo "--- muscle"; tail -c 120 logs/muscle.stderr.log; echo
echo "--- mafft"; tail -c 120 logs/mafft.stderr.log; echo
ls -la aln
ps -eo pid,etime,pcpu,rss,cmd --sort=-pcpu | head -6 | cut -c1-110
free -m | head -2
date -u
