# session cell 159 | 2026-10-18 12:51 | env: r | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/sOTU_Bd_panelA/data/panelA_isolate_scores.csv","sha256":"fe9f135fd29dc4dc8ac9833488f869dfb5a9bc6e6b8dddc0cb24aca4070b25a9"}]

cd amphib_R; mkdir -p ../sOTU_Bd_panelA/data ../sOTU_Bd_panelA/output
Rscript -e '
source("R/00_functions.R")
X  <- fread(file.path(DER, "amphibac_isolate_matches.csv")); st <- fread(file.path(DER, "sOTU_stats.csv"))
E <- X[Bd_function %in% c("Inhibitory","Non-inhibitory","Facilitating") & sOTUs_3mm != "", .(sOTU = unlist(strsplit(sOTUs_3mm, ";"))), by = .(FastaID, Bd_function, inhib_score, Reference, host_norm)]
Es <- E[!is.na(inhib_score)]; keep <- Es[, .N, by=sOTU][N>=5, sOTU]
D <- merge(Es[sOTU %in% keep], st[, .(sOTU, mean_rel_pct, prevalence_pct)], by="sOTU")
ids <- D[, .(m=median(inhib_score)), by=sOTU][order(m)][, sOTU_id := sprintf("sOTU_%03d", .I)]
D <- merge(D, ids[, .(sOTU, sOTU_id)], by="sOTU")
setnames(D, c("sOTU","host_norm"), c("sOTU_V4_sequence","host_species"))
setcolorder(D, c("sOTU_id","sOTU_V4_sequence","FastaID","Reference","host_species","Bd_function","inhib_score","mean_rel_pct","prevalence_pct"))
fwrite(D[order(sOTU_id, FastaID)], "../sOTU_Bd_panelA/data/panelA_isolate_scores.csv")
cat(nrow(D), uniqueN(D$sOTU_id), uniqueN(D$FastaID), "\n")'
head -c 600 ../sOTU_Bd_panelA/data/panelA_isolate_scores.csv