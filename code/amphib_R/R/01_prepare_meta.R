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
M <- as(as.matrix(dt), "CsparseMatrix"); rm(dt); invisible(gc())
rownames(M) <- ids
totals <- colSums(M)
fwrite(data.table(SampleID = names(totals), total_reads = totals, included = totals >= MIN_READS),
       file.path(DER, "sample_totals.csv"))
keep <- totals >= MIN_READS
R <- t(t(M[, keep]) / totals[keep]) * 100             # relative abundance (%)
R <- drop0(R); R <- R[rowSums(R) > 0, ]
st <- data.table(sOTU = rownames(R), reads = rowSums(M[rownames(R), keep]))
nz <- tabulate(R@i + 1, nbins = nrow(R))
st[, n_samples := nz]
st[, prevalence_pct := 100 * n_samples / sum(keep)]
st[, mean_rel_pct := rowMeans(R)]
mx <- numeric(nrow(R)); tm <- tapply(R@x, R@i + 1, max); mx[as.integer(names(tm))] <- tm
st[, max_rel_pct := mx]
st[, mean_when_present_pct := mean_rel_pct / (prevalence_pct / 100)]
st[, rank := frank(-mean_rel_pct, ties.method = "min")]
fwrite(st, file.path(DER, "sOTU_stats.csv"))
saveRDS(R, file.path(DER, "rel_abund_sparse.rds"))
message("samples used: ", sum(keep), " of ", length(keep), "; sOTUs present: ", nrow(R))
invisible(file.remove(txt))
