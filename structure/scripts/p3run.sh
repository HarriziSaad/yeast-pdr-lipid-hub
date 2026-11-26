#!/bin/bash
export MAMBA_ROOT_PREFIX=/home/pdr/micromamba
mkdir -p /home/pdr/work/phase3; cd /home/pdr/work/phase3
exec /home/pdr/bin/micromamba run -n pdr "$@"
