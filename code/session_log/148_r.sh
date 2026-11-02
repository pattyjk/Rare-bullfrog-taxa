# session cell 148 | 2026-10-18 12:08 | env: r | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/figures/rank_abundance.png","sha256":"7b6e6f74e7e3803f5248629df6eeb62435fa43aecbc2f4b59617e60045442bb6","preview_sha256":"a67dd3fe04886fcfc8c5448274719b6d23b471362e90db19860c13ca011c68c4","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/figures/cultured_vs_abundance.png","sha256":"4ced7010f1d96ce91796426ccbf806ce74d8abcae50cee4a839c7a7c02a184b7","preview_sha256":"5c53f05f54ad336e5c8133b5b879d3a6c8923032c8690c5473b87dcdca213978","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/detection_vs_function.csv","sha256":"160631d0128244b6e844fa735e1f7f52ba73424f18f26da62a79899ab1fdb082"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/function_vs_abundance_tests.csv","sha256":"e8492539b09f4ce68e1325289abfc2f5ab4d8e484eb7a8fbbb56e38bc8c600a7"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/per_sample_cultured_share.csv","sha256":"c22eb24f73efdfa4057ee813dbea85a61e5cc9336cbbb2729b510e46083528cb"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/near_match_thresholds.csv","sha256":"7ccaf845b9909610dc24048b81e877ad7c141f17b09d331fde7159792c991686"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/cultured_by_abundance_class.csv","sha256":"c35d27d192de1a2a57e77568fad0e3c566408022fb764cd510be650024…[+249 chars]

cd amphib_R
cat > R/06_cultured_share.R <<'EOF'
## 06. How abundant are cultured (AmphiBac) taxa in the amphibian skin microbiome? ----
## Inputs : derived/sOTU_stats.csv, rel_abund_sparse.rds (01); derived/amphibac_isolate_matches.csv (04);
##          output/sOTU_Bd_function_variance.csv (05)
## Outputs: output/near_match_thresholds.csv, cultured_by_abundance_class.csv, per_sample_cultured_share.csv,
##          detection_vs_function.csv, function_vs_abundance_tests.csv;
##          figures/rank_abundance.png, figures/cultured_vs_abundance.png
source("R/00_functions.R")
st <- fread(file.path(DER, "sOTU_stats.csv")); R <- readRDS(file.path(DER, "rel_abund_sparse.rds"))
X  <- fread(file.path(DER, "amphibac_isolate_matches.csv")); S <- fread(file.path(OUT, "sOTU_Bd_function_variance.csv"))
st[, abundance_class := abund_class(mean_rel_pct)]

## near-match thresholds ------------------------------------------------------------
v4n <- X[has_V4 == TRUE | detection_exact == "detected", .N]
thr <- rbindlist(lapply(c(0, 1, 2, 3, 5), function(t) {
  det <- X$detection_exact == "detected" | (!is.na(X$nearest_mismatch) & X$nearest_mismatch <= t)
  data.table(max_mismatches = t, identity = sprintf(">=%.1f%%", 100 * (90 - t) / 90), isolates_detected = sum(det), of_V4_spanning = v4n,
             pct_detected = 100 * sum(det) / v4n)
}))
fwrite(thr, file.path(OUT, "near_match_thresholds.csv")); print(thr, digits = 3)

## cultured share by abundance class ---------------------------------------------------
cult_exact <- unique(unlist(strsplit(X[exact_sOTUs != "", exact_sOTUs], ";")))
cult_3 <- unique(c(cult_exact, X[sOTUs_3mm != "", sOTUs_3mm]))
cf <- st[, .(n_sOTUs = .N, n_cultured_exact = sum(sOTU %in% cult_exact), n_cultured_3mm = sum(sOTU %in% cult_3)), by = abundance_class][order(abundance_class)]
cf[, `:=`(pct_cultured_exact = 100 * n_cultured_exact / n_sOTUs, pct_cultured_3mm = 100 * n_cultured_3mm / n_sOTUs)]
fwrite(cf, file.path(OUT, "cultured_by_abundance_class.csv")); print(cf, digits = 3)

## per-sample shares (<=3 mismatches) -------------------------------------------------
P <- R > 0; inC <- rownames(R) %in% cult_3
ps <- data.table(SampleID = colnames(R), richness = colSums(P), cultured_sOTUs = colSums(P[inC, ]),
                 pct_reads_cultured = colSums(R[inC, ]))[, pct_sOTUs_cultured := 100 * cultured_sOTUs / richness]
fwrite(ps, file.path(OUT, "per_sample_cultured_share.csv"))
message("samples with zero cultured sOTUs: ", ps[cultured_sOTUs == 0, .N], " | median % sOTUs: ", round(median(ps$pct_sOTUs_cultured), 1),
        " | median % reads: ", round(median(ps$pct_reads_cultured), 1), " | Spearman n vs %reads: ",
        round(cor(ps$cultured_sOTUs, ps$pct_reads_cultured, method = "spearman"), 2))

## Bd function vs abundance (sOTU level, <=3 mm) ----------------------------------------
inh <- S[function_class == "Inhibitory"]; non <- S[function_class == "No effect / mixed"]
Xb <- X[!is.na(Bd_function) & detection_3mm != "unknown (no V4 coverage)"]
dtab <- table(Xb$detection_3mm, Xb$Bd_function == "Inhibitory")
tests <- data.table(test = c("Mann-Whitney mean rel. abundance: inhibitory vs no-effect/mixed sOTUs",
                             "Mann-Whitney prevalence: inhibitory vs no-effect/mixed sOTUs",
                             "Spearman: fraction inhibitory vs mean rel. abundance (sOTUs)",
                             "Chi-square: isolate inhibitory vs detected/not detected (<=3 mm)"),
                    stat = c(wilcox.test(inh$mean_rel_pct, non$mean_rel_pct)$statistic, wilcox.test(inh$prevalence_pct, non$prevalence_pct)$statistic,
                             cor(S$frac_inhib, S$mean_rel_pct, method = "spearman"), chisq.test(dtab)$statistic),
                    p = c(wilcox.test(inh$mean_rel_pct, non$mean_rel_pct)$p.value, wilcox.test(inh$prevalence_pct, non$prevalence_pct)$p.value,
                          cor.test(S$frac_inhib, S$mean_rel_pct, method = "spearman", exact = FALSE)$p.value, chisq.test(dtab)$p.value))
fwrite(tests, file.path(OUT, "function_vs_abundance_tests.csv"))
fwrite(as.data.table(table(detection = Xb$detection_3mm, Bd_function = Xb$Bd_function)), file.path(OUT, "detection_vs_function.csv"))
print(tests, digits = 3)

## figures -------------------------------------------------------------------------------
st2 <- merge(st, S[, .(sOTU, function_class)], by = "sOTU", all.x = TRUE)
pa <- ggplot(st2, aes(rank, mean_rel_pct)) + geom_line(colour = "grey70", linewidth = 0.6) +
  geom_point(data = st2[!is.na(function_class)][order(function_class == "Inhibitory")], aes(colour = function_class), size = 0.8, alpha = 0.8) +
  geom_hline(yintercept = 0.1, linetype = 3, linewidth = 0.3) + scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) + scale_colour_manual(values = COL_FUN) +
  labs(x = "sOTU abundance rank", y = "Mean relative abundance (%)", title = "Cultured isolates mostly match top-ranked sOTUs") +
  theme_pub() + theme(legend.position = c(0.25, 0.2))
pci <- S[, .(n = .N, k = sum(function_class == "Inhibitory")), by = abundance_class][order(abundance_class)][, pct := 100 * k / n]
pb <- ggplot(pci, aes(abundance_class, pct)) + geom_col(fill = COL_FUN["Inhibitory"], width = 0.65) + geom_text(aes(label = paste0(k, "/", n)), vjust = -0.4, size = 2) +
  scale_x_discrete(labels = c("Dominant\n\u22651%", "Interm.\n0.1\u20131%", "Rare\n<0.1%")) + ylim(0, 75) +
  labs(x = NULL, y = "Bd-inhibitory sOTUs (%)", title = "Rare sOTUs not more often inhibitory") + theme_pub()
ggsave(file.path(FIG, "rank_abundance.png"), pa + pb + plot_layout(widths = c(2, 1)) + plot_annotation(tag_levels = "a"), width = 7.2, height = 3.2, dpi = 300)
cl <- melt(cf, id.vars = c("abundance_class","n_sOTUs"), measure.vars = c("pct_cultured_exact","pct_cultured_3mm"), variable.name = "match", value.name = "pct")
cl[, match := factor(match, labels = c("Exact V4 match", "\u22643 mismatches (\u226596.7%)"))]
p1 <- ggplot(cl, aes(abundance_class, pct, fill = match)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  scale_fill_manual(values = c("#2166AC", "#92C5DE")) + scale_x_discrete(labels = c("Dominant\n\u22651%", "Interm.\n0.1\u20131%", "Rare\n<0.1%")) +
  labs(x = NULL, y = "sOTUs with a cultured\nAmphiBac match (%)", title = "All dominant sOTUs are cultured") + theme_pub() + theme(legend.position = "top")
psl <- melt(ps, id.vars = "SampleID", measure.vars = c("pct_sOTUs_cultured","pct_reads_cultured"))[, variable := factor(variable, labels = c("sOTUs", "Reads"))]
p2 <- ggplot(psl, aes(variable, value)) + geom_jitter(width = 0.18, size = 0.2, alpha = 0.15) + geom_boxplot(outlier.shape = NA, width = 0.5, fill = NA) +
  labs(x = NULL, y = "Share matched to cultured\nisolates per sample (%)", title = "Few sOTUs, most reads") + theme_pub()
bins <- ps[, .(med = median(pct_reads_cultured), .N), by = .(bin = floor(cultured_sOTUs / 10) * 10 + 5)][N >= 20][order(bin)]
p3 <- ggplot(ps, aes(cultured_sOTUs, pct_reads_cultured)) + geom_point(size = 0.3, alpha = 0.25, colour = "#2166AC") +
  geom_line(data = bins, aes(bin, med)) + geom_point(data = bins, aes(bin, med), size = 0.8) +
  labs(x = "Cultured sOTUs per sample (n)", y = "Reads matched to cultured\nisolates (%)", title = "Read share declines above ~20 matches") + theme_pub()
ggsave(file.path(FIG, "cultured_vs_abundance.png"), p1 + p2 + p3 + plot_layout(widths = c(1.1, 0.8, 1.1)) + plot_annotation(tag_levels = "a"), width = 7.4, height = 3, dpi = 300)
EOF
time Rscript R/06_cultured_share.R 2>&1 | grep -v "^Warning\|^ *$" | tail -30