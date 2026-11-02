# session cell 145 | 2026-10-18 12:06 | env: r | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/figures/sOTU_Bd_variance.png","sha256":"68999747eb85cc80f46534bcf74465cd4d58f047fa211e982eb72351f425a93c","preview_sha256":"a20734ba8b55e7edfe1ce54f6b95c2b4f996166317b8aadcd1f52dc13497d41a","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/sOTU_Bd_function_variance.csv","sha256":"4e0cf702cf9b535808f31112957bc3daff679135a95687c6e584cdda4167e5be"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/Bd_variance_partitioning.csv","sha256":"2394ad8f158ca9c31e18e97d2e683953003bdf0da287f79b33034d2a37443681"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/R/05_sOTU_function_variance.R","sha256":"5c902bb24adacfdfc91d60feacb194a94e7f9af3d85078e27c64748898102cfe"}]

cd amphib_R
cat > R/05_sOTU_function_variance.R <<'EOF'
## 05. Per-sOTU Bd function and variability of Bd inhibition among isolates ----
## Inputs : derived/amphibac_isolate_matches.csv (04), derived/sOTU_stats.csv (01)
## For every meta-analysis sOTU matched (<=3 mm) by >=1 Bd-tested isolate:
##   function counts, fraction inhibitory, label agreement (% isolates sharing the majority label),
##   inhibition-score n, mean, median, SD, IQR, MAD, min, max, range; number of studies contributing.
## Variance partitioning: one-way ANOVA ICC of inhibition score with sOTU as group (sOTUs with >=2 scored isolates)
##   = share of isolate-level score variance that lies between sOTUs; same for study (Reference) within sOTU.
## Outputs: output/sOTU_Bd_function_variance.csv, output/Bd_variance_partitioning.csv, figures/sOTU_Bd_variance.png
source("R/00_functions.R")
X  <- fread(file.path(DER, "amphibac_isolate_matches.csv"))
st <- fread(file.path(DER, "sOTU_stats.csv"))
E <- X[!is.na(Bd_function) & sOTUs_3mm != "", .(sOTU = unlist(strsplit(sOTUs_3mm, ";"))), by = .(FastaID, Bd_function, inhib_score, Reference, host_norm)]
lab_major <- function(f) { t <- table(f); names(t)[which.max(t)] }
S <- E[, .(n_isolates = uniqueN(FastaID), n_inhib = sum(Bd_function == "Inhibitory"), n_non = sum(Bd_function == "Non-inhibitory"),
           n_fac = sum(Bd_function == "Facilitating"), majority_label = lab_major(Bd_function),
           label_agreement_pct = 100 * max(table(Bd_function)) / .N, n_studies = uniqueN(Reference), n_hosts = uniqueN(host_norm),
           n_scored = sum(!is.na(inhib_score)), mean_score = mean(inhib_score, na.rm = TRUE), median_score = median(inhib_score, na.rm = TRUE),
           sd_score = sd(inhib_score, na.rm = TRUE), iqr_score = IQR(inhib_score, na.rm = TRUE), mad_score = mad(inhib_score, na.rm = TRUE),
           min_score = suppressWarnings(min(inhib_score, na.rm = TRUE)), max_score = suppressWarnings(max(inhib_score, na.rm = TRUE))), by = sOTU]
S[n_scored == 0, c("mean_score","median_score","min_score","max_score") := NA]
S[, range_score := max_score - min_score][, frac_inhib := n_inhib / n_isolates]
S[, function_class := fifelse(frac_inhib >= 0.5, "Inhibitory", fifelse(n_fac / n_isolates >= 0.5, "Enhancing", "No effect / mixed"))]
S <- merge(S, st[, .(sOTU, mean_rel_pct, prevalence_pct, mean_when_present_pct)], by = "sOTU")
S[, abundance_class := abund_class(mean_rel_pct)]
fwrite(S[order(-n_isolates)], file.path(OUT, "sOTU_Bd_function_variance.csv"))

icc <- function(y, g) { ok <- !is.na(y); y <- y[ok]; g <- factor(g[ok]); keep <- g %in% names(which(table(g) >= 2)); y <- y[keep]; g <- droplevels(g[keep])
  a <- anova(lm(y ~ g)); k <- nlevels(g); N <- length(y); n0 <- (N - sum(table(g)^2) / N) / (k - 1)
  msb <- a$`Mean Sq`[1]; msw <- a$`Mean Sq`[2]; c(groups = k, isolates = N, ICC = (msb - msw) / (msb + (n0 - 1) * msw), F_p = a$`Pr(>F)`[1]) }
Es <- E[!is.na(inhib_score)]
vp <- rbind(sOTU = icc(Es$inhib_score, Es$sOTU), study = icc(Es$inhib_score, Es$Reference),
            sOTU_x_study = icc(Es$inhib_score, paste(Es$sOTU, Es$Reference)))
# within-sOTU: does study explain score differences among isolates of the same sOTU? (sOTUs with >=2 studies)
multi <- Es[, if (uniqueN(Reference) >= 2 & .N >= 5) .SD, by = sOTU]
within_study_p <- multi[, .(p = tryCatch(kruskal.test(inhib_score, factor(Reference))$p.value, error = function(e) NA_real_)), by = sOTU]
vp <- data.table(grouping = rownames(vp), vp)
vp <- rbind(vp, data.table(grouping = "study within sOTU: % of multi-study sOTUs with Kruskal p<0.05", groups = nrow(within_study_p),
                           isolates = nrow(multi), ICC = 100 * mean(within_study_p$p < 0.05, na.rm = TRUE), F_p = NA))
fwrite(vp, file.path(OUT, "Bd_variance_partitioning.csv"))
print(vp, digits = 3)

## figure -------------------------------------------------------------------------
V <- S[n_scored >= 5][order(median_score)][, idx := .I]
q <- Es[sOTU %in% V$sOTU, .(q25 = quantile(inhib_score, 0.25), q75 = quantile(inhib_score, 0.75)), by = sOTU]
V <- merge(V, q, by = "sOTU")
colc <- c("Dominant (>=1%)" = "#08306B", "Intermediate (0.1-1%)" = "#4292C6", "Rare (<0.1%)" = "#BDBDBD")
pa <- ggplot(V, aes(idx, median_score, colour = abundance_class)) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey40") +
  geom_linerange(aes(ymin = pmax(min_score, -250), ymax = pmin(max_score, 200)), linewidth = 0.2, alpha = 0.35) +
  geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 0.6) + geom_point(size = 0.6, colour = "black") +
  scale_colour_manual(values = colc) + coord_cartesian(ylim = c(-250, 200)) +
  labs(x = sprintf("sOTUs with >=5 scored isolates (n = %d), ranked by median", nrow(V)), y = "Bd inhibition score (%)",
       title = "Many sOTUs span inhibition to enhancement") + theme_pub() + theme(legend.position = c(0.25, 0.85), axis.text.x = element_blank(), axis.ticks.x = element_blank())
pb <- ggplot(V, aes(abundance_class, sd_score)) + geom_boxplot(outlier.shape = NA, fill = "#D1E5F0", width = 0.6) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.5, colour = "#2166AC") +
  scale_x_discrete(labels = c("Dominant\n>=1%", "Interm.\n0.1-1%", "Rare\n<0.1%")) +
  labs(x = NULL, y = "SD of inhibition score\namong isolates (%)", title = "Spread is similar across classes") + theme_pub()
pc <- ggplot(S[n_isolates >= 3], aes(label_agreement_pct)) + geom_histogram(binwidth = 10, boundary = 0, fill = "#4D4D4D", colour = "white") +
  labs(x = "Isolates sharing the majority\nBd label (%)", y = "sOTUs (n)", title = "Label agreement within sOTUs") + theme_pub()
p <- pa + pb + pc + plot_layout(widths = c(2, 1, 1)) + plot_annotation(tag_levels = "a")
ggsave(file.path(FIG, "sOTU_Bd_variance.png"), p, width = 7.4, height = 2.9, dpi = 300)
message("sOTUs with >=5 scored isolates: ", nrow(V), " | median SD: ", round(median(V$sd_score), 1),
        " | % spanning inhibition (>50) and enhancement (<0): ", round(100 * mean(V$max_score > 50 & V$min_score < 0), 1))
print(V[, .(n = .N, median_sd = median(sd_score), median_iqr = median(iqr_score)), by = abundance_class][order(abundance_class)], digits = 3)
print(kruskal.test(sd_score ~ abundance_class, data = V))
print(S[n_isolates >= 3, .(n = .N, median_agreement = median(label_agreement_pct), pct_unanimous = 100 * mean(label_agreement_pct == 100))])
EOF
time Rscript R/05_sOTU_function_variance.R 2>&1 | tail -30