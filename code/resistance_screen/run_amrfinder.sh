#!/usr/bin/env bash
# AMRFinderPlus (--plus) on BF1-BF7 Prokka outputs.
# Usage: bash run_amrfinder.sh <fg_inputs_dir with BFn/BFn.{faa,fna,gff}> <outdir>
# Needs: conda/bioconda ncbi-amrfinderplus (tested 4.2.7) and blastp (bundled with it).
set -euo pipefail
IN=$1; OUT=$2; mkdir -p "$OUT"
DB="$OUT/amrdb"; [ -d "$DB/latest" ] || amrfinder_update -d "$DB"      # database 2026-10-29.1 was used
for i in 1 2 3 4 5 6 7; do
  d=$IN/BF$i; faa=$d/BF$i.faa
  if [ $i = 3 ]; then
    # IHHBNBAB_01063 crashes amr_report (HMM-domain parser, AMRFinderPlus 4.2.7); run without it,
    # then call it by BLASTP against the AMRFinderPlus reference proteins.
    awk '/^>/{skip=($1==">IHHBNBAB_01063")} !skip' $faa > "$OUT/BF3_minus_IHHBNBAB_01063.faa"
    awk '/^>/{keep=($1==">IHHBNBAB_01063")} keep' $faa > "$OUT/IHHBNBAB_01063.faa"
    blastp -query "$OUT/IHHBNBAB_01063.faa" -subject "$DB/latest/AMRProt.fa" -evalue 1e-10 \
      -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" -max_target_seqs 5 > "$OUT/IHHBNBAB_01063_blastp.tsv"
    faa="$OUT/BF3_minus_IHHBNBAB_01063.faa"
  fi
  amrfinder -p $faa -n $d/BF$i.fna -g $d/BF$i.gff -a prokka --plus -d "$DB/latest" --threads 4 \
            --name BF$i -o "$OUT/BF$i.amrfinder.tsv" > "$OUT/BF$i.amrfinder.log" 2>&1
done
amrfinder --database_version -d "$DB/latest" > "$OUT/amrfinder_versions.txt" 2>&1 || true
