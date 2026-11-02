# session cell 153 | 2026-10-18 12:36 | env: r | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/run_all.R","sha256":"1ed57040f97758ec026e31c78b57ed9bf8e79d9e81150aea022494e9bd71df80"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/blast_bullfrog_asvs.sh","sha256":"75f0b662b227aa7229d24106fccf4d26f50ec8bec0ebcfee9581b463e08e8c59"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/R/optional_fetch_ncbi_refs.R","sha256":"b7ed5a7e099ac5f0181bd1ba28deb6cf9c6e70b3c44b0d0e9754dcf7d5ca0288"}]

cd amphib_R
cat > R/optional_fetch_ncbi_refs.R <<'EOF'
## OPTIONAL. Re-download the NCBI 16S reference sequences used for exact-V4 matching -------
## The copies used in the analysis (downloaded 2026-10-16) are already in data/raw/ncbi_refs/.
## Re-running gives a newer GenBank snapshot, so counts may differ slightly. Requires internet access.
## Query per genus: "<Genus>"[Organism] AND (16S ribosomal RNA[Title] OR 16S rRNA[Title]) AND 1000:1700[SLEN]; first 1,500 records.
source("R/00_functions.R")
fetch_genus <- function(g, outdir = file.path(RAW, "ncbi_refs"), retmax = 1500) {
  base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  term <- utils::URLencode(sprintf('"%s"[Organism] AND (16S ribosomal RNA[Title] OR 16S rRNA[Title]) AND 1000:1700[SLEN]', g), reserved = TRUE)
  x <- paste(readLines(sprintf("%sesearch.fcgi?db=nuccore&term=%s&retmax=%d&usehistory=y", base, term, retmax), warn = FALSE), collapse = "")
  n <- min(as.integer(sub(".*?<Count>(\\d+)</Count>.*", "\\1", x)), retmax)
  we <- sub(".*<WebEnv>(.*?)</WebEnv>.*", "\\1", x); qk <- sub(".*<QueryKey>(\\d+)</QueryKey>.*", "\\1", x)
  fa <- unlist(lapply(seq(0, n - 1, by = 500), function(s) { Sys.sleep(0.4)
    readLines(sprintf("%sefetch.fcgi?db=nuccore&query_key=%s&WebEnv=%s&rettype=fasta&retmode=text&retstart=%d&retmax=500", base, qk, we, s), warn = FALSE) }))
  writeLines(fa, file.path(outdir, paste0(g, ".fasta"))); message(g, ": ", sum(grepl("^>", fa)), " sequences")
}
if (sys.nframe() == 0) for (g in REF_GENERA) fetch_genus(g)
EOF
cat > blast_bullfrog_asvs.sh <<'EOF'
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
EOF
chmod +x blast_bullfrog_asvs.sh
cat > run_all.R <<'EOF'
## Run the full pipeline from the project root: Rscript run_all.R
## Steps 01 and 04 are the slow ones (~1 min each); their outputs are already in data/derived/, so you can
## start from any later script. Set RERUN_PREP <- FALSE to skip 01 and 04 and reuse data/derived/.
RERUN_PREP <- TRUE
steps <- c("R/01_prepare_meta.R", "R/02_meta_genus_abundance.R", "R/03_bullfrog_genus_abundance.R", "R/04_amphibac_match.R",
           "R/05_sOTU_function_variance.R", "R/06_cultured_share.R", "R/07_abundance_class_function.R", "R/08_source_host.R")
if (!RERUN_PREP) steps <- setdiff(steps, c("R/01_prepare_meta.R", "R/04_amphibac_match.R"))
for (s in steps) { message("\n===== ", s, " ====="); t0 <- Sys.time(); source(s, local = new.env()); message("done in ", round(difftime(Sys.time(), t0, units = "secs")), " s") }
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
EOF
echo written