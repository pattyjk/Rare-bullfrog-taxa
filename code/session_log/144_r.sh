# session cell 144 | 2026-10-18 11:50 | env: r | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/data/derived/amphibac_isolate_matches.csv","sha256":"ce2d808de4ee84d1782187dc8487c04fae64f56e97f7c82717dbf5274185e9d7"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/R/04_amphibac_match.R","sha256":"a3fdc362eff2b6b46517cbf04f6026868884831b8ca0b5296ef2af30777d9bf9"}]

cd amphib_R
cat > R/04_amphibac_match.R <<'EOF'
## 04. Match AmphiBac isolates to amphibian meta-analysis sOTUs ----------------
## Inputs : AmphibBac_Metadata_2025.1-tested.tsv, AmphibBac_2025.1-tested.fasta (Woodhams et al. 2015 database, 2025.1),
##          derived/sOTU_stats.csv (01)
## Matching: (1) exact - a 90-bp sOTU with reads in the amphibian data occurs exactly within the isolate 16S (either strand);
##           (2) near - isolates spanning V4 (90 bp after 515F) but without an exact match are assigned the nearest sOTU
##               if <=3 mismatches (>=96.7% identity; allows for Sanger error). Ties -> most abundant sOTU.
##           Only sOTUs with reads are used: the taxonomy file also contains primer-containing / reverse-oriented
##           sequences with zero reads that would give spurious matches.
## Bd function: Bd.Inhibition normalised (Inhibitory / Non-inhibitory / Facilitating). Inhibition score =
##           (1 - Proportional.Growth.Bd) x 100 (96-well CFS assays only; 100 = full inhibition, <0 = enhanced growth).
## Output : data/derived/amphibac_isolate_matches.csv
source("R/00_functions.R")
st <- fread(file.path(DER, "sOTU_stats.csv"))
pres <- st$sOTU
meta <- fread(file.path(RAW, "AmphibBac_Metadata_2025.1-tested.tsv"), sep = "\t", quote = "", colClasses = "character")
fa <- read_fasta(file.path(RAW, "AmphibBac_2025.1-tested.fasta"))
kmers <- function(s, k = 90) { n <- nchar(s); if (n < k) character(0) else substring(s, 1:(n - k + 1), k:n) }
pres_set <- pres
exact <- vapply(seq_along(fa), function(i) {
  s <- fa[[i]]; h <- unique(c(kmers(s), kmers(revcomp(s))))
  paste(sort(h[h %in% pres_set]), collapse = ";")
}, "")
names(exact) <- names(fa)
win <- v4_window(unname(fa), 90); names(win) <- names(fa)
need <- names(fa)[exact == "" & !is.na(win) & grepl("^[ACGTN]+$", win)]
PM <- seq_matrix(pres); ab <- st$mean_rel_pct
near <- rbindlist(lapply(need, function(k) {
  d <- hamming_to(win[[k]], PM); dm <- min(d); ii <- which(d == dm)
  data.table(FastaID = k, nearest_mismatch = dm, n_tied = length(ii), nearest_sOTU = pres[ii[which.max(ab[ii])]])
}))
X <- data.table(FastaID = names(fa), exact_sOTUs = exact, has_V4 = !is.na(win))
X <- merge(X, near, by = "FastaID", all.x = TRUE)
X[, sOTUs_3mm := fifelse(exact_sOTUs != "", exact_sOTUs, fifelse(!is.na(nearest_mismatch) & nearest_mismatch <= 3, nearest_sOTU, ""))]
X[, detection_exact := fifelse(exact_sOTUs != "", "detected", fifelse(has_V4, "not detected (V4 not in amphibian data)", "unknown (no V4 coverage)"))]
X[, detection_3mm := fifelse(sOTUs_3mm != "", "detected", fifelse(has_V4, "not detected (V4 not in amphibian data)", "unknown (no V4 coverage)"))]
fixg <- c(gephryomantis = "gephyromantis", ischconemia = "ischnocnema", hypsiboas = "boana")
norm_host <- function(x) { w <- strsplit(tolower(trimws(gsub("[._]", " ", x))), "\\s+")
  vapply(w, function(z) { if (!length(z)) return(NA_character_); z[1] <- ifelse(z[1] %in% names(fixg), fixg[z[1]], z[1]); paste(head(z, 2), collapse = " ") }, "") }
m <- meta[, .(FastaID, Isolate.ID, Host.Species, host_norm = norm_host(Host.Species), Country, Life.Stage, Reference, Assay.Method,
              Bd.Inhibition, Bd_function = norm_bd(Bd.Inhibition), Proportional.Growth.Bd)]
m[, PG := suppressWarnings(as.numeric(Proportional.Growth.Bd))][, inhib_score := (1 - PG) * 100]
X <- merge(m, X, by = "FastaID")
fwrite(X, file.path(DER, "amphibac_isolate_matches.csv"))
message("isolates: ", nrow(X), " | detected exact: ", X[detection_exact == "detected", .N], " | detected <=3 mm: ", X[detection_3mm == "detected", .N],
        " | V4 spanning: ", X[has_V4 == TRUE, .N])
print(table(X$nearest_mismatch)[1:6])
EOF
time Rscript R/04_amphibac_match.R 2>&1 | tail -6