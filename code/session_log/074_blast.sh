# session cell 74 | 2026-10-17 16:57 | env: blast | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/ref16S_11genera_ncbi.tar.gz","sha256":"c6ca11279e95fff9e73dc516a539af1e4e8e39ce0b9f1d8d2e37397b491ef470"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/sal/asv_blast.tsv","sha256":"8709f3e28864b5cf4b0c9edabdd836a8757000f09bedae8f4f9d67eb30934566"}]

tar czf ref16S_11genera_ncbi.tar.gz -C meta ref; ls -la ref16S_11genera_ncbi.tar.gz | awk '{print $5}'
cd meta/db; export BLASTDB=$(pwd); blastdbcmd -db 16S_ribosomal_RNA -info | grep -E "sequences"
blastn -task blastn -query ../../sal/asvs.fasta -db 16S_ribosomal_RNA -max_target_seqs 50 -num_threads 10 -evalue 1e-10 \
 -outfmt "6 qseqid sacc pident length mismatch gapopen qstart qend evalue bitscore staxids sscinames stitle" > ../../sal/asv_blast.tsv
cut -f1 ../../sal/asv_blast.tsv | sort -u | wc -l