#!/usr/bin/env bash
echo "# run date (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"
python -c "import sys,pyhmmer,pandas;print('python',sys.version.split()[0],'| pyhmmer',pyhmmer.__version__,'| pandas',pandas.__version__)" 2>&1 || true
(hmmsearch -h 2>&1 | sed -n 2p) || true
(cd-hit -h 2>&1 | grep -m1 "CD-HIT version") || true
(snakemake --version 2>&1 | sed 's/^/snakemake /') || true
echo "# authoritative versions: project/tools/software_versions.tsv"
exit 0
