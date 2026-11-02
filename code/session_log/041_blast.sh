# session cell 41 | 2026-10-15 21:01 | env: blast | language: bash
# files written: []

cd meta/db; export BLASTDB=$(pwd); blastdbcmd -db 16S_ribosomal_RNA -entry all -outfmt "%S" | grep -iE "^(Psychromicrobium|Arthrobacter|Cryobacterium|Ancylobacter|Rhodococcoides)" | awk '{print $1}' | sort | uniq -c