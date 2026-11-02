## 03. Target-genus relative abundance in the bullfrog ASV table --------------
## Inputs : bullfrog_asv_table.txt (42 samples x 3,700 ASVs), bullfrog_dna-sequences.fasta (120-bp V4 after 515F),
##          bullfrog_asv_blast_16SRefSeq.tsv (blastn 2.17.0 -task blastn vs NCBI 16S_ribosomal_RNA, evalue 1e-10,
##          max_target_seqs 50; see README), derived/ncbi_ref_v4.csv (from 02), meta_analysis_taxonomy.txt
## Methods: BLAST97_top (top hit >=97%, aln >=100 bp); BLAST97_tied (genus among hits tied at best bitscore);
##          exact120 (identical to 120 bp after 515F of an NCBI sequence); exact120_min2refs; exact120_unique_min2refs;
##          metaTaxonomy90_97 (first 90 bp identical to a meta-analysis sOTU assigned to the genus at >=97%);
##          any_method = union of BLAST97_top, BLAST97_tied, exact120, metaTaxonomy90_97.
## Outputs: output/bullfrog_genus_summary.csv, bullfrog_genus_relabund_by_sample.csv, bullfrog_genus_relabund_long.csv,
##          bullfrog_target_ASVs.csv
source("R/00_functions.R")
tab <- fread(file.path(RAW, "bullfrog_asv_table.txt")); setnames(tab, 1, "ASV")
cnt <- as.matrix(tab[, -1]); rownames(cnt) <- tab$ASV
seqs <- read_fasta(file.path(RAW, "bullfrog_dna-sequences.fasta"))[rownames(cnt)]
tot <- colSums(cnt); inc <- tot >= MIN_READS
rel <- t(t(cnt) / tot) * 100
message("samples: ", ncol(cnt), " (", sum(!inc), " with <", MIN_READS, " reads excluded from summaries)")
lab <- function(g) fifelse(g == "Rhodococcus", "Rhodococcus (other)", g)
genus_of <- function(sci) {
  w <- tstrsplit(sci, " ", keep = 1:2)
  fifelse(w[[1]] %in% c("Rhodococcus","Rhodococcoides") & w[[2]] %in% RHODOCOCCOIDES_EPITHETS, "Rhodococcoides",
  fifelse(w[[1]] == "Rhodococcus", "Rhodococcus (other)", w[[1]]))
}
b <- fread(file.path(RAW, "bullfrog_asv_blast_16SRefSeq.tsv"), sep = "\t", quote = "",
           col.names = c("qseqid","sacc","pident","length","mismatch","gapopen","qstart","qend","evalue","bitscore","staxids","sscinames","stitle"))
b[, genus := genus_of(sscinames)]
top  <- b[, .SD[1], by = qseqid]
tied <- b[, .SD[bitscore == max(bitscore)], by = qseqid][, .(tied_genera = paste(sort(unique(genus)), collapse = ";"), n_tied = .N), by = qseqid]
A <- data.table(ASV = rownames(cnt), seq = seqs)
A <- merge(A, top[, .(ASV = qseqid, blast_top_hit = sscinames, blast_genus = genus, pident, aln = length)], by = "ASV", all.x = TRUE)
A <- merge(A, tied[, .(ASV = qseqid, tied_genera, n_tied)], by = "ASV", all.x = TRUE)
refs <- fread(file.path(DER, "ncbi_ref_v4.csv"))[!is.na(v4_120)]
r120 <- refs[, .(n_refs = .N), by = .(v4_120, tg = lab(ref_genus))]
r120all <- refs[, .(exact120_genera = paste(sort(unique(lab(ref_genus))), collapse = ";"), exact120_n_refs = .N,
                    exact120_species = paste(head(names(sort(table(species), decreasing = TRUE)), 4), collapse = "; ")), by = v4_120]
A <- merge(A, r120all, by.x = "seq", by.y = "v4_120", all.x = TRUE)
tax <- unique(fread(file.path(RAW, "meta_analysis_taxonomy.txt"), sep = "\t", quote = ""), by = "sOTU")
tax[, genus_tax := sub(".*; ", "", sub(";\\s*$", "", Taxonomy))][, c("hg","he") := tstrsplit(Hit_name, " ", keep = 1:2)]
tax[, meta_genus := fifelse(hg == "Rhodococcus" & he %in% RHODOCOCCOIDES_EPITHETS, "Rhodococcoides",
                    fifelse(genus_tax == "Rhodococcus", "Rhodococcus (other)", genus_tax))]
A[, p90 := substr(seq, 1, 90)]
A <- merge(A, tax[Per_ident >= 97, .(p90 = sOTU, meta_genus, meta_hit = Hit_name)], by = "p90", all.x = TRUE)
has <- function(x, g) vapply(strsplit(fifelse(is.na(x), "", x), ";"), function(z) g %in% z, TRUE)
M <- list()
for (g in TARGETS) {
  m <- list(
    BLAST97_top  = !is.na(A$blast_genus) & A$blast_genus == g & A$pident >= 97 & A$aln >= 100,
    BLAST97_tied = has(A$tied_genera, g) & !is.na(A$pident) & A$pident >= 97,
    exact120     = A$seq %in% r120[tg == g, v4_120],
    exact120_min2refs = A$seq %in% r120[tg == g & n_refs >= 2, v4_120],
    metaTaxonomy90_97 = !is.na(A$meta_genus) & A$meta_genus == g)
  m$exact120_unique_min2refs <- m$exact120_min2refs & !is.na(A$exact120_genera) & A$exact120_genera == g
  m$any_method <- m$BLAST97_top | m$BLAST97_tied | m$exact120 | m$metaTaxonomy90_97
  M[[g]] <- m
}
cols <- list(); summ <- list()
for (g in TARGETS) for (k in names(M[[g]])) {
  ids <- A$ASV[M[[g]][[k]]]
  r <- if (length(ids)) colSums(rel[ids, , drop = FALSE]) else setNames(rep(0, ncol(rel)), colnames(rel))
  nm <- paste0(sub(" \\(other\\)", "_other", g), "__", k); cols[[nm]] <- r; ri <- r[inc]
  summ[[nm]] <- data.table(genus = g, method = k, n_ASVs = length(ids), prevalence_pct = 100 * mean(ri > 0), n_samples = sum(ri > 0),
                           mean_rel_pct = mean(ri), median_rel_pct = median(ri), max_rel_pct = max(ri))
}
summ <- rbindlist(summ); fwrite(summ, file.path(OUT, "bullfrog_genus_summary.csv"))
W <- data.table(SampleID = colnames(rel), total_reads = tot, included_in_summary = inc, as.data.table(cols))
fwrite(W, file.path(OUT, "bullfrog_genus_relabund_by_sample.csv"))
L <- melt(W, id.vars = c("SampleID","total_reads","included_in_summary"), variable.name = "gm", value.name = "rel_abund_pct")[, c("genus","method") := tstrsplit(gm, "__")][, gm := NULL]
fwrite(L, file.path(OUT, "bullfrog_genus_relabund_long.csv"))
flag <- rbindlist(lapply(TARGETS, function(g) data.table(ASV = A$ASV, genus = g, methods = apply(sapply(M[[g]], identity), 1, function(z) paste(names(z)[z & names(z) != "any_method"], collapse = ";")))[methods != ""]))
flag <- merge(flag, A[, .(ASV, blast_top_hit, pident, tied_genera, exact120_genera, exact120_n_refs, exact120_species, meta_hit, seq)], by = "ASV")
flag[, `:=`(mean_rel_pct = rowMeans(rel[ASV, inc, drop = FALSE]), prevalence_pct = 100 * rowMeans(rel[ASV, inc, drop = FALSE] > 0))]
fwrite(flag[order(genus, -mean_rel_pct)], file.path(OUT, "bullfrog_target_ASVs.csv"))
print(dcast(summ, genus ~ method, value.var = "mean_rel_pct"), digits = 3)
