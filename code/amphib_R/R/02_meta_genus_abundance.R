## 02. Target-genus relative abundance in the amphibian meta-analysis ----------
## Inputs : meta_analysis_taxonomy.txt (BLAST top hit to NCBI 16S RefSeq for each sOTU),
##          ncbi_refs/*.fasta (NCBI nuccore 16S, 1000-1700 bp, per genus; downloaded 2026-10-16),
##          derived rel_abund_sparse.rds, amphib_only_map.txt
## Methods: BLAST97 / BLASTany = genus of top hit (>=97% / any identity);
##          exactV4 = sOTU identical to the 90 bp after 515F in >=1 NCBI sequence of the genus;
##          exactV4_min2refs = supported by >=2 NCBI records (single records are often misidentified);
##          exactV4_unique_min2refs = additionally matches no other genus in the 11-genus reference set.
## Outputs: derived/ncbi_ref_v4.csv; output/meta_genus_summary.csv, meta_genus_by_host_metadata.csv,
##          meta_genus_relabund_by_sample.csv, meta_target_sOTUs.csv
source("R/00_functions.R")
R  <- readRDS(file.path(DER, "rel_abund_sparse.rds"))
mp <- fread(file.path(RAW, "amphib_only_map.txt"))[`#SampleID` %in% colnames(R)]
setnames(mp, "#SampleID", "SampleID")

## reference V4 table (also used by 03) ----------------------------------------
refs <- rbindlist(lapply(REF_GENERA, function(g) {
  s <- read_fasta(file.path(RAW, "ncbi_refs", paste0(g, ".fasta")))
  w <- strsplit(names(s), " ")
  data.table(ref_genus = g, accession = vapply(w, `[`, "", 1),
             species = vapply(w, function(z) paste(z[2:3], collapse = " "), ""),
             v4_90 = v4_window(s, 90), v4_120 = v4_window(s, 120))
}))
refs[!grepl("^[ACGT]+$", v4_90), v4_90 := NA][!grepl("^[ACGT]+$", v4_120), v4_120 := NA]
fwrite(refs, file.path(DER, "ncbi_ref_v4.csv"))
message("reference sequences with V4: ", refs[!is.na(v4_90), .N])
lab <- function(g) fifelse(g == "Rhodococcus", "Rhodococcus (other)", g)

## BLAST-based assignment from the meta-analysis taxonomy file ------------------
tax <- unique(fread(file.path(RAW, "meta_analysis_taxonomy.txt"), sep = "\t", quote = ""), by = "sOTU")
tax[, genus_tax := sub(".*; ", "", sub(";\\s*$", "", Taxonomy))]
tax[, c("hit_genus","hit_epithet") := tstrsplit(Hit_name, " ", keep = 1:2)]
tax[, target := NA_character_]
tax[hit_genus %in% c("Rhodococcus","Rhodococcoides") & hit_epithet %in% RHODOCOCCOIDES_EPITHETS, target := "Rhodococcoides"]
for (g in setdiff(TARGETS, c("Rhodococcoides","Rhodococcus (other)")))
  tax[is.na(target) & (grepl(paste0("; ", g, "(;|$)"), Taxonomy) | hit_genus == g), target := g]
tax[is.na(target) & grepl("; Rhodococcus(;|$)", Taxonomy), target := "Rhodococcus (other)"]

## exact V4 matches -------------------------------------------------------------
ex <- refs[!is.na(v4_90), .(n_refs = .N), by = .(v4_90, tg = lab(ref_genus))]
ex_all <- refs[!is.na(v4_90), .(all_genera = paste(sort(unique(lab(ref_genus))), collapse = ";")), by = v4_90]
ex <- merge(ex, ex_all, by = "v4_90")[v4_90 %in% rownames(R)]

sets <- list()
for (g in TARGETS) {
  sets[[paste0(g, "__BLAST97")]]  <- intersect(tax[target == g & Per_ident >= 97, sOTU], rownames(R))
  sets[[paste0(g, "__BLASTany")]] <- intersect(tax[target == g, sOTU], rownames(R))
  if (g != "Rhodococcus (other)") {
    e <- ex[tg == g]
    sets[[paste0(g, "__exactV4")]] <- e$v4_90
    sets[[paste0(g, "__exactV4_min2refs")]] <- e[n_refs >= 2, v4_90]
    sets[[paste0(g, "__exactV4_unique_min2refs")]] <- e[n_refs >= 2 & all_genera == g, v4_90]
  }
}
G <- sapply(sets, function(s) if (length(s)) colSums(R[s, , drop = FALSE]) else rep(0, ncol(R)))
summ <- rbindlist(lapply(names(sets), function(k) {
  r <- G[, k]; p <- r[r > 0]
  data.table(genus = sub("__.*", "", k), method = sub(".*__", "", k), n_sOTUs = length(sets[[k]]),
             prevalence_pct = 100 * mean(r > 0), n_samples = sum(r > 0), mean_rel_pct = mean(r),
             median_rel_when_present_pct = if (length(p)) median(p) else NA, max_rel_pct = max(r))
}))
fwrite(summ, file.path(OUT, "meta_genus_summary.csv"))
wide <- data.table(SampleID = colnames(R), G)
fwrite(merge(mp[, .(SampleID, host_order, host_family, host_genus, host_species, lifestage, country, biogeo_realm, latitude_deg, longitude_deg, doi)], wide, by = "SampleID"),
       file.path(OUT, "meta_genus_relabund_by_sample.csv"))
strata <- c("host_order","host_family","lifestage","biogeo_realm","country")
long <- melt(wide, id.vars = "SampleID", variable.name = "gm", value.name = "rel")[, c("genus","method") := tstrsplit(gm, "__")]
long <- merge(long, mp[, c("SampleID", strata), with = FALSE], by = "SampleID")
bystr <- rbindlist(lapply(strata, function(s) long[, .(n_samples = .N, prevalence_pct = 100 * mean(rel > 0), mean_rel_pct = mean(rel)),
                                                  by = c(s, "genus", "method")][n_samples >= 20][, stratum := s][, level := get(s)][, (s) := NULL]))
fwrite(bystr, file.path(OUT, "meta_genus_by_host_metadata.csv"))
st <- fread(file.path(DER, "sOTU_stats.csv"))
det <- rbindlist(lapply(names(sets), function(k) if (length(sets[[k]])) data.table(sOTU = sets[[k]], genus = sub("__.*", "", k), method = sub(".*__", "", k))))
det <- det[, .(methods = paste(method, collapse = ";")), by = .(sOTU, genus)]
det <- merge(det, tax[, .(sOTU, blast_top_hit = Hit_name, blast_pct_identity = Per_ident)], by = "sOTU", all.x = TRUE)
det <- merge(det, ex_all[, .(sOTU = v4_90, exactV4_ref_genera = all_genera)], by = "sOTU", all.x = TRUE)
det <- merge(det, st[, .(sOTU, prevalence_pct, mean_rel_pct, max_rel_pct)], by = "sOTU")
fwrite(det[order(genus, -mean_rel_pct)], file.path(OUT, "meta_target_sOTUs.csv"))
print(dcast(summ, genus ~ method, value.var = "prevalence_pct"), digits = 3)
