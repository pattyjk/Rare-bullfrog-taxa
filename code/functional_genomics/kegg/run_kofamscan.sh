#!/usr/bin/env bash
# KofamScan on BF1-BF7 Prokka proteomes (conda env 'annot'). CPU capped at 3 (GNU parallel -j 3, hmmsearch --cpu 1).
set -euo pipefail
DB=${DB:-..}   # directory holding profiles/ and ko_list (NOT distributed; download below)
# curl -O https://www.genome.jp/ftp/db/kofam/profiles.tar.gz && tar xzf profiles.tar.gz   # Last-Modified 2026-10-27
# curl -O https://www.genome.jp/ftp/db/kofam/ko_list.gz && gunzip ko_list.gz              # Last-Modified 2026-10-26
for i in BF1 BF2 BF3 BF4 BF5 BF6 BF7; do sed "s/^>/>${i}__/" fg_inputs/$i/$i.faa; done > all7.faa
mkdir -p kofam_tmp
exec_annotation -p $DB/profiles -k $DB/ko_list --cpu 3 --tmp-dir kofam_tmp -f detail-tsv -o kofam_all7_detail.tsv all7.faa
# Module definitions: fetched from https://rest.kegg.jp/list/module and /get/<M> (batches of 10) on 2026-11-02 -> modules_raw/
python scripts/kegg_postprocess.py
