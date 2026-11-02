# session cell 39 | 2026-10-15 21:01 | env: blast | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/asvs.fasta","sha256":"8b002585a02ad97086cbba30d747f003dfd8b41ddbeb9890c3cfa009eda6be83"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/asv_blast.tsv","sha256":"aa0fe06dded6d160ef999429c5ffab49b9652d58fa626a6a754127f255fa239a"}]

cd meta/db; cp "<ARTIFACT dna-sequences.fasta>" ../asvs.fasta
export BLASTDB=$(pwd)
blastn -version | head -1
blastdbcmd -db 16S_ribosomal_RNA -info | grep -E "sequences|Date"
blastn -task blastn -query ../asvs.fasta -db 16S_ribosomal_RNA -max_target_seqs 50 -num_threads 10 -evalue 1e-10 \
 -outfmt "6 qseqid sacc pident length mismatch gapopen qstart qend evalue bitscore staxids sscinames stitle" > ../asv_blast.tsv
wc -l ../asv_blast.tsv; cut -f1 ../asv_blast.tsv | sort -u | wc -l; head -2 ../asv_blast.tsv