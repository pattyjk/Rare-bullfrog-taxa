# session cell 152 | 2026-10-18 12:35 | env: r | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/figures/abundance_class_definitions.png","sha256":"af8e5c62a6754fd05af27cd4edd766b6ed3b8ba8ce3e05d6d6e85e39ef1b80a6","preview_sha256":"699795bd5523d5f2af90724977c2ad2a4b194738f52052aa8581ecb31acefa6a","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/figures/abundance_class_function_isolates.png","sha256":"9781db78149a05e9677e4ffe7a4fb120a14d2fdd2a34f24d72828bc6997b8ad8","preview_sha256":"9bea74ef6f84a1f7499cde2dd987432c1abe034f83d9b42390a1b0f13e8410b8","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/figures/source_host_abundance.png","sha256":"0a5c1b21dc1a06fff8378f5973404b98297b25e5c1e33d960fa2381ab3586866","preview_sha256":"18adfddb4589472a7696cd5bc70bdfe8731eb79cf2136fc5e1d2baf32bbd1e4e","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/figures/abundance_class_function_taxa.png","sha256":"691a1571eba83edb138e94f5211acd664d2cc1e1b1402a2caab380fd59249d11","preview_sha256":"c02fbdd50a82939246f39fb146157f44e4fd46178c672851414d45fd18d9320c","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/source_host_summary.csv","sha256":"89d8cad3d2b2327138f6c191ea5c4fc086d06dd4620fae884586dd15db350fa8"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/abundance_class_tests.csv","sha256":"5eaa77d45ca9cbe682447db557861b5326feba8f91cee576c23e…[+1513 chars]

cd amphib_R
sed -i.bak 's/labs(x = NULL, y = "Units by Bd function (%)") + theme_pub() + theme(legend.position = "none")/labs(x = NULL, y = "Units by Bd function (%)") + theme_pub() + theme(legend.position = "top")/' R/07_abundance_class_function.R && rm R/*.bak
Rscript R/07_abundance_class_function.R >/dev/null 2>&1
cat > R/08_source_host.R <<'EOF'
## 08. Are cultured taxa more abundant on the host species they were isolated from? ----
## Inputs : derived/amphibac_isolate_matches.csv (04), derived/rel_abund_sparse.rds (01), amphib_only_map.txt
## Source host = AmphiBac Host.Species; if not sampled in the meta-analysis, closest relative: same host genus, else same
##   family (families for genera absent from the meta-analysis assigned manually below). Unit = matched sOTU set x source host.
## Comparisons: (1) pooled - mean in source-host samples vs mean over all other samples;
##              (2) host-level - source-host group mean vs each other host group (species/genus/family) with >=5 samples;
##                  pct_other_groups_below = % of other groups with strictly lower mean.
## Outputs: output/source_host_units.csv, source_host_summary.csv; figures/source_host_abundance.png
source("R/00_functions.R")
X <- fread(file.path(DER, "amphibac_isolate_matches.csv"))[sOTUs_3mm != ""]
R <- readRDS(file.path(DER, "rel_abund_sparse.rds"))
mp <- fread(file.path(RAW, "amphib_only_map.txt"))[`#SampleID` %in% colnames(R)]; setnames(mp, "#SampleID", "SampleID")
mp <- mp[match(colnames(R), SampleID)][, `:=`(hg = tolower(trimws(host_genus)), hf = tolower(trimws(host_family)))][, hsp := paste(hg, tolower(trimws(host_species)))]
g2f <- mp[, .N, by = .(hg, hf)][order(-N)][!duplicated(hg), setNames(hf, hg)]
fam_manual <- c(salamandra="salamandridae", gephyromantis="mantellidae", boana="hylidae", brachycephalus="brachycephalidae", litoria="hylidae",
  espadarana="centrolenidae", silverstoneia="dendrobatidae", hyalinobatrachium="centrolenidae", anodontohyla="microhylidae", hemidactylium="plethodontidae",
  pristimantis="craugastoridae", stumffia="microhylidae", ischnocnema="brachycephalidae", bufo="bufonidae", haddadus="craugastoridae", dendrobates="dendrobatidae",
  leiopelma="leiopelmatidae", plethodonthyla="microhylidae", allobates="aromobatidae", hylodes="hylodidae", paradoxohyla="microhylidae", telmatobius="telmatobiidae",
  psychrophrynella="craugastoridae", adenomera="leptodactylidae", bryophryne="craugastoridae", gastrotheca="hemiphractidae", oreobates="craugastoridae", osteocephalus="hylidae")
X[, hg := sub(" .*", "", host_norm)]
X[, hf := fifelse(hg %in% names(g2f), g2f[hg], fam_manual[hg])]
X[, home_level := fifelse(host_norm %in% mp$hsp & !grepl(" (sp|aff|cf)$", host_norm), "same host species",
                  fifelse(hg %in% mp$hg, "same host genus", fifelse(hf %in% mp$hf, "same host family", "no relative in meta-analysis")))]
print(X[, .N, by = home_level])
X[, fn := fifelse(Bd_function == "Inhibitory", "Inhibitory", fifelse(Bd_function %in% c("Non-inhibitory","Facilitating"), "Non-inhibitory/facilitating", NA_character_))]
U <- X[home_level != "no relative in meta-analysis", .(n_isolates = .N, n_bd = sum(!is.na(fn)), n_inhib = sum(fn == "Inhibitory", na.rm = TRUE),
       hg = hg[1], hf = hf[1], home_level = home_level[1]), by = .(sOTUs_3mm, host_norm)]
U[, function_class := fifelse(n_bd == 0, "untested", fifelse(n_inhib / n_bd >= 0.5, "Inhibitory", "Non-inhibitory/facilitating"))]
res <- rbindlist(lapply(seq_len(nrow(U)), function(i) {
  u <- U[i]; ab <- colSums(R[strsplit(u$sOTUs_3mm, ";")[[1]], , drop = FALSE])
  key <- switch(u$home_level, "same host species" = "hsp", "same host genus" = "hg", "same host family" = "hf")
  grp <- mp[[key]]; me <- switch(key, hsp = u$host_norm, hg = u$hg, hf = u$hf); home <- grp == me
  gm <- tapply(ab, grp, mean); gn <- table(grp); gm <- gm[names(gn)[gn >= 5]]
  oth <- gm[names(gm) != me]; hgm <- if (me %in% names(gm)) gm[[me]] else NA
  data.table(home_n = sum(home), home_mean = mean(ab[home]), other_mean = mean(ab[!home]), home_prev = 100 * mean(ab[home] > 0), other_prev = 100 * mean(ab[!home] > 0),
             home_group_mean = hgm, median_other_groups = median(oth), pct_other_groups_below = if (is.na(hgm)) NA else 100 * mean(oth < hgm),
             home_group_rank = if (is.na(hgm)) NA else rank(-gm, ties.method = "min")[[me]], n_groups = length(gm))
}))
U <- cbind(U, res)
fwrite(U, file.path(OUT, "source_host_units.csv"))
summ <- U[, .(n_units = .N, pooled_median_source = median(home_mean), pooled_median_other = median(other_mean), pooled_pct_source_higher = 100 * mean(home_mean > other_mean),
              n_hostlevel = sum(!is.na(home_group_mean)), pct_source_above_median_other = 100 * mean(home_group_mean > median_other_groups, na.rm = TRUE),
              median_pct_groups_below = median(pct_other_groups_below, na.rm = TRUE), pct_absent_on_source = 100 * mean(home_group_mean == 0, na.rm = TRUE),
              source_is_top = sum(home_group_rank == 1, na.rm = TRUE)), by = home_level]
fwrite(summ, file.path(OUT, "source_host_summary.csv")); print(summ, digits = 3)
sp <- U[home_level == "same host species" & !is.na(pct_other_groups_below) & function_class != "untested"]
print(sp[, .(n = .N, median_pct_below = median(pct_other_groups_below), pct_absent = 100 * mean(home_group_mean == 0)), by = function_class])
print(wilcox.test(pct_other_groups_below ~ function_class, data = sp[home_group_mean > 0]))
pc <- 1e-4; f <- function(v) log10(v + pc); tk <- c(0, 1e-3, 1e-2, 1e-1, 1, 10)
D <- U[!is.na(home_group_mean)][, home_level := factor(home_level, levels = c("same host family","same host genus","same host species"))]
pa <- ggplot(D, aes(f(home_group_mean), f(median_other_groups), colour = home_level)) + geom_abline(linetype = 2, linewidth = 0.3) + geom_point(size = 0.5, alpha = 0.5) +
  scale_colour_manual(values = c("#BDBDBD", "#67A9CF", "#2166AC")) + scale_x_continuous(breaks = f(tk), labels = tk) + scale_y_continuous(breaks = f(tk), labels = tk) +
  labs(x = "Mean relative abundance on source host (%)", y = "Median across other host groups (%)", title = "Source host vs other hosts") + theme_pub() + theme(legend.position = c(0.25, 0.88))
pb <- ggplot(sp, aes(pct_other_groups_below, colour = function_class)) + geom_freqpoly(binwidth = 10, boundary = 0, linewidth = 0.8) +
  scale_colour_manual(values = c("Inhibitory" = "#2166AC", "Non-inhibitory/facilitating" = "#E08214")) +
  labs(x = "Other host species with lower abundance\nthan the source host (%)", y = "Isolate-host units (n)", title = "Source host: often top-ranked or taxon absent") +
  theme_pub() + theme(legend.position = "top")
ggsave(file.path(FIG, "source_host_abundance.png"), pa + pb + plot_annotation(tag_levels = "a"), width = 7, height = 3.3, dpi = 300)
EOF
time Rscript R/08_source_host.R 2>&1 | grep -v "^Warning\|Removed" | tail -25