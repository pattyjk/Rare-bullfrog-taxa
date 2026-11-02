#!/usr/bin/env bash
## Re-create data/raw/bullfrog_asv_blast_16SRefSeq.tsv (already provided). Requires BLAST+ (tested with 2.17.0).
## Database: NCBI 16S_ribosomal_RNA (27,647 sequences, 2026-10-15 release) + taxdb for scientific names.
set -euo pipefail
mkdir -p blastdb && cd blastdb
curl -O https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz && tar xzf 16S_ribosomal_RNA.tar.gz
curl -O https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && tar xzf taxdb.tar.gz
export BLASTDB=$(pwd); cd ..
blastn -task blastn -query data/raw/bullfrog_dna-sequences.fasta -db 16S_ribosomal_RNA -max_target_seqs 50 -evalue 1e-10 -num_threads 8 \
  -outfmt "6 qseqid sacc pident length mismatch gapopen qstart qend evalue bitscore staxids sscinames stitle" \
  > data/raw/bullfrog_asv_blast_16SRefSeq.tsv
