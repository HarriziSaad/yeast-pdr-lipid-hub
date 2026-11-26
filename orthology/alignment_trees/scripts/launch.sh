#!/bin/bash
log=$1; pidf=$2; shift 2
setsid nohup "$@" > "$log" 2>&1 < /dev/null &
echo $! > "$pidf"
echo "launched pid $(cat $pidf): $*"
