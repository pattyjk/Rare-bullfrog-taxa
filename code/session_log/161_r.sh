# session cell 161 | 2026-10-18 12:52 | env: r | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/sOTU_Bd_panelA/panelA_sOTU_Bd_spread.R","sha256":"2bf46cea0cf46ff0abdfe92779c873a2a59a0be6fd09a4e98c11e29a55249326"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/sOTU_Bd_panelA/output/panelA_sOTU_Bd_spread.png","sha256":"bf5eaf428707c36fa54e44ce4f54e143c414bf461bbf866a6ceafba58539dd9c","preview_sha256":"bf5eaf428707c36fa54e44ce4f54e143c414bf461bbf866a6ceafba58539dd9c","preview_content_type":"image/png"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/sOTU_Bd_panelA/output/panelA_sOTU_Bd_spread.pdf","sha256":"7cf893dd7dcbd480ee2cd532bad6de2c7461d60a2ce9cbe9089dcae59f603945"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/sOTU_Bd_panelA/output/panelA_sOTU_summary.csv","sha256":"2a62f4503b3cb60590b5867e339a911886fed605c8eb16fd2b865de47232e713"}]

cd sOTU_Bd_panelA
cat > panelA_sOTU_Bd_spread.R <<'EOF'
## Within-sOTU spread of Bd inhibition among AmphiBac isolates (stand-alone)
## Run from this folder: Rscript panelA_sOTU_Bd_spread.R
## Input : data/panelA_isolate_scores.csv  (one row per scored isolate; see README/methods)
## Output: output/panelA_sOTU_summary.csv, output/panelA_sOTU_Bd_spread.png/.pdf
## Packages: data.table, ggplot2
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
dir.create("output", showWarnings = FALSE)
D <- fread("data/panelA_isolate_scores.csv")

abund_class <- function(p) factor(fifelse(p >= 1, "Dominant (>=1%)", fifelse(p >= 0.1, "Intermediate (0.1-1%)", "Rare (<0.1%)")),
                                  levels = c("Dominant (>=1%)", "Intermediate (0.1-1%)", "Rare (<0.1%)"))
S <- D[, .(n_scored = .N, n_studies = uniqueN(Reference), n_hosts = uniqueN(host_species),
           n_inhib = sum(Bd_function == "Inhibitory"), n_non = sum(Bd_function == "Non-inhibitory"), n_fac = sum(Bd_function == "Facilitating"),
           median_score = median(inhib_score), mean_score = mean(inhib_score), sd_score = sd(inhib_score),
           q25 = quantile(inhib_score, 0.25), q75 = quantile(inhib_score, 0.75),
           min_score = min(inhib_score), max_score = max(inhib_score), mean_rel_pct = mean_rel_pct[1], prevalence_pct = prevalence_pct[1]),
       by = .(sOTU_id, sOTU_V4_sequence)]
S <- S[order(median_score)][, idx := .I][, abundance_class := abund_class(mean_rel_pct)]
S[, spans_inhib_and_enh := max_score > 50 & min_score < 0]
fwrite(S, "output/panelA_sOTU_summary.csv")

colc <- c("Dominant (>=1%)" = "#08306B", "Intermediate (0.1-1%)" = "#4292C6", "Rare (<0.1%)" = "#BDBDBD")
clab <- c("Dominant (\u22651%)", "Intermediate (0.1\u20131%)", "Rare (<0.1%)")
YL <- c(-250, 200)   # display window; whiskers beyond it are clipped (true min/max are in the summary CSV)
p <- ggplot(S, aes(idx, median_score, colour = abundance_class)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey40") +
  geom_linerange(aes(ymin = pmax(min_score, YL[1]), ymax = pmin(max_score, YL[2])), linewidth = 0.25, alpha = 0.35) +
  geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 0.7) +
  geom_point(size = 0.7, colour = "black") +
  scale_colour_manual(values = colc, labels = clab, name = "Mean relative abundance") +
  coord_cartesian(ylim = YL) +
  labs(x = sprintf("sOTUs with \u22655 scored isolates (n = %d), ranked by median", nrow(S)), y = "Bd inhibition score (%)") +
  theme_classic(base_size = 9) +
  theme(legend.position = c(0.25, 0.85), legend.title = element_text(size = 8), legend.text = element_text(size = 7.5),
        legend.key.height = unit(0.35, "cm"), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("output/panelA_sOTU_Bd_spread.png", p, width = 4.2, height = 3, dpi = 300)
ggsave("output/panelA_sOTU_Bd_spread.pdf", p, width = 4.2, height = 3)

message(sprintf("sOTUs: %d | isolates: %d | studies: %d | median within-sOTU SD: %.1f | IQR of SDs: %.1f-%.1f | spanning >50 and <0: %.1f%% | median score <0: %d",
  nrow(S), nrow(D), uniqueN(D$Reference), median(S$sd_score), quantile(S$sd_score, .25), quantile(S$sd_score, .75),
  100 * mean(S$spans_inhib_and_enh), sum(S$median_score < 0)))
print(S[, .(n = .N, median_sd = round(median(sd_score), 1)), by = abundance_class][order(abundance_class)])
EOF
Rscript panelA_sOTU_Bd_spread.R 2>&1; Rscript -e 'library(data.table);d<-fread("data/panelA_isolate_scores.csv");print(table(d$Bd_function));print(range(d$inhib_score));s<-fread("output/panelA_sOTU_summary.csv");cat(sum(s$min_score< -250), sum(s$max_score>200), median(s$n_scored), max(s$n_scored), "\n")'