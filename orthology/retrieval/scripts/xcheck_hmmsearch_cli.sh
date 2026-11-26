#!/usr/bin/env bash
set -u
faa=$1; hmmdir=$2; hits=$3; out=$4; thr=${5:-10}
mkdir -p "$out"; : > "$out/cli_domains.tsv"
for acc in PF00005 PF01061 PF06422 PF14510 PF19055; do
  hmmsearch --cut_ga --cpu $thr --noali --domtblout "$out/$acc.domtbl" -o /dev/null "$hmmdir/$acc.hmm" "$faa"
  grep -v '^#' "$out/$acc.domtbl" | awk -v a=$acc '{print $1"\t"a"\t"$18"\t"$19}' >> "$out/cli_domains.tsv"
done
sort "$out/cli_domains.tsv" > "$out/cli_sorted.tsv"
awk -F'\t' 'NR>1{print $1"\t"$5"\t"$11"\t"$12}' "$hits" | sort > "$out/pyhmmer_sorted.tsv"
echo "CLI domains: $(wc -l < $out/cli_sorted.tsv)  pyhmmer domains: $(wc -l < $out/pyhmmer_sorted.tsv)  identical rows (seq, family, ali_from, ali_to): $(comm -12 $out/cli_sorted.tsv $out/pyhmmer_sorted.tsv | wc -l)  only CLI: $(comm -23 $out/cli_sorted.tsv $out/pyhmmer_sorted.tsv | wc -l)  only pyhmmer: $(comm -13 $out/cli_sorted.tsv $out/pyhmmer_sorted.tsv | wc -l)"
