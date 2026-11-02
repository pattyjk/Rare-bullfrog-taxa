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
