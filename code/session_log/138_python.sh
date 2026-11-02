# session cell 138 | 2026-10-18 11:46 | env: python | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/R/00_functions.R","sha256":"21c6f14e3b151ea1535333ecac2d7ca3ec0f072cc36c2478fff39a40b3ad5bd1"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/R/01_prepare_meta.R","sha256":"3fb438adfed98379e8af64e160ac17d2115d5fd913978c67660f76e806a1e2a1"}]

cd amphib_R
cat > R/00_functions.R <<'EOF'
## Shared paths, constants and helper functions --------------------------------
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(ggplot2); library(patchwork)
})
RAW <- "data/raw"; DER <- "data/derived"; FIG <- "figures"; OUT <- "output"
for (d in c(DER, FIG, OUT)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

TARGETS <- c("Rhodococcoides","Psychromicrobium","Arthrobacter","Cryobacterium",
             "Ancylobacter","Sphingomonas","Kineococcus","Rhodococcus (other)")
REF_GENERA <- c("Rhodococcoides","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter",
                "Sphingomonas","Kineococcus","Rhodococcus","Paenarthrobacter","Pseudarthrobacter","Glutamicibacter")
# species moved from Rhodococcus to Rhodococcoides (Val-Calvo & Vazquez-Boland 2024), old and new epithets
RHODOCOCCOIDES_EPITHETS <- c("fascians","corynebacterioides","kroppenstedtii","kyotonensis","kyotonense",
                             "navarretei","trifolii","yunnanensis","yunnanense","luteus")
PRIMER_515F <- "GTG[CT]CAGC[AC]GCCGCGGTAA"
PRIMER_806R_RC <- "ATTAGA[AT]ACCC[CGT][ACGT]GTAGTCC"
MIN_READS <- 1000
COL_FUN <- c("Inhibitory"="#2166AC","No effect / mixed"="#E08214","Enhancing"="#7B3294")
CLASS_LEVELS <- c("Dominant (>=1%)","Intermediate (0.1-1%)","Rare (<0.1%)")

read_fasta <- function(path) {
  x <- readLines(path, warn = FALSE)
  h <- grepl("^>", x); id <- cumsum(h)
  seqs <- vapply(split(x[!h], id[!h]), function(s) toupper(paste(s, collapse = "")), "")
  nm <- sub("^>", "", x[h]); names(seqs) <- nm[as.integer(names(seqs))]
  gsub("-", "", seqs)
}
revcomp <- function(s) chartr("ACGTRYKMBDHVN", "TGCAYRMKVHDBN", stringi_rev(s))
stringi_rev <- function(s) vapply(strsplit(s, NULL), function(z) paste(rev(z), collapse = ""), "")
# V4 window: n bp immediately after the 515F primer (either strand); NA if not covered
v4_window <- function(s, n = 90) {
  out <- rep(NA_character_, length(s))
  for (strand in 1:2) {
    q <- if (strand == 1) s else revcomp(s)
    m <- regexpr(PRIMER_515F, q, perl = TRUE)
    ok <- is.na(out) & m > 0 & (m + attr(m, "match.length") - 1 + n) <= nchar(q)
    out[ok] <- substr(q[ok], m[ok] + attr(m, "match.length")[ok], m[ok] + attr(m, "match.length")[ok] + n - 1)
  }
  out
}
abund_class <- function(x) factor(ifelse(x >= 1, CLASS_LEVELS[1], ifelse(x >= 0.1, CLASS_LEVELS[2], CLASS_LEVELS[3])), levels = CLASS_LEVELS)
# Hamming distance of one sequence to every row of a character matrix (one base per column)
hamming_to <- function(q, M) rowSums(M != matrix(strsplit(q, NULL)[[1]], nrow(M), ncol(M), byrow = TRUE))
seq_matrix <- function(s) do.call(rbind, strsplit(s, NULL))
norm_bd <- function(x) {
  x <- trimws(x)
  fifelse(x %in% c("Inhibitory","inhibitory","Inhibitor"), "Inhibitory",
  fifelse(x %in% c("NonInhibitory","NoEffect","No Effect","No effect","Non-inhibitor","ns"), "Non-inhibitory",
  fifelse(x %in% c("Facilitating","enhancing"), "Facilitating", NA_character_)))
}
theme_pub <- function() theme_classic(base_size = 8) +
  theme(plot.title = element_text(size = 8, hjust = 0), axis.text = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 6), legend.title = element_blank(), legend.key.size = unit(3, "mm"),
        plot.tag = element_text(face = "bold", size = 9))
EOF

cat > R/01_prepare_meta.R <<'EOF'
## 01. Parse the amphibian meta-analysis OTU table -----------------------------
## Input : data/raw/amphib_only_otu_table.zip (sOTU x sample counts; 90-bp V4 deblur sOTUs as row IDs)
## Output: data/derived/sample_totals.csv, sOTU_stats.csv, rel_abund_sparse.rds
## Notes : 52 malformed rows (wrong field count) are dropped with awk; numeric-ID rows (all zero) are dropped.
source("R/00_functions.R")
zipf <- file.path(RAW, "amphib_only_otu_table.zip")
txt <- file.path(DER, "amphib_only_otu_table.txt")
if (!file.exists(txt)) unzip(zipf, exdir = DER)
hdr <- strsplit(readLines(txt, n = 1), "\t")[[1]]          # header has no first-column name
nf <- length(hdr) + 1
cmd <- sprintf("awk -F'\\t' 'NR>1 && NF==%d && length($1)==90' %s", nf, shQuote(txt))
dt <- fread(cmd = cmd, header = FALSE, sep = "\t")
setnames(dt, c("sOTU", hdr))
message("sOTU rows kept: ", nrow(dt))
ids <- dt$sOTU; dt[, sOTU := NULL]
M <- as(as.matrix(dt), "CsparseMatrix"); rm(dt); gc()
rownames(M) <- ids
totals <- colSums(M)
fwrite(data.table(SampleID = names(totals), total_reads = totals, included = totals >= MIN_READS),
       file.path(DER, "sample_totals.csv"))
keep <- totals >= MIN_READS
R <- t(t(M[, keep]) / totals[keep]) * 100             # relative abundance (%)
R <- drop0(R); R <- R[rowSums(R) > 0, ]
st <- data.table(sOTU = rownames(R), reads = rowSums(M[rownames(R), keep]),
                 n_samples = diff(R@p)[0] , stringsAsFactors = FALSE)
nz <- tabulate(R@i + 1, nbins = nrow(R))
st[, n_samples := nz]
st[, prevalence_pct := 100 * n_samples / sum(keep)]
st[, mean_rel_pct := rowMeans(R)]
st[, max_rel_pct := apply(R, 1, max)]
st[, mean_when_present_pct := mean_rel_pct / (prevalence_pct / 100)]
st[, rank := frank(-mean_rel_pct, ties.method = "min")]
fwrite(st, file.path(DER, "sOTU_stats.csv"))
saveRDS(R, file.path(DER, "rel_abund_sparse.rds"))
message("samples used: ", sum(keep), " of ", length(keep), "; sOTUs present: ", nrow(R))
file.remove(txt)
EOF
echo ok