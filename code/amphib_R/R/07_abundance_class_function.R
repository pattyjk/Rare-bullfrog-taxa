## 07. Bd function and inhibition strength by abundance class (three class definitions) ----
## Inputs : derived/amphibac_isolate_matches.csv (04), output/sOTU_Bd_function_variance.csv (05),
##          derived/rel_abund_sparse.rds (01), amphib_only_map.txt
## Class definitions (thresholds dominant >=1%, intermediate 0.1-1%, rare <0.1% in all cases):
##   global       = sOTU mean relative abundance across all 2,411 samples (zeros included)
##   when_present = mean across samples where the sOTU is detected
##   host         = mean on the isolate's own source host species (host species with >=5 samples); unit = sOTU x host;
##                  'Not detected' = zero reads in all samples of that host
## Function per unit: Inhibitory if >=50% of Bd-tested isolates inhibit, Enhancing if >=50% facilitate, else No effect/mixed.
## Outputs: output/abundance_class_summary.csv, abundance_class_tests.csv, host_specific_units.csv;
##          figures/abundance_class_function_taxa.png, abundance_class_function_isolates.png, abundance_class_definitions.png
source("R/00_functions.R")
X <- fread(file.path(DER, "amphibac_isolate_matches.csv"))[Bd_function %in% c("Inhibitory","Non-inhibitory","Facilitating") & sOTUs_3mm != ""]
S <- fread(file.path(OUT, "sOTU_Bd_function_variance.csv"))
R <- readRDS(file.path(DER, "rel_abund_sparse.rds"))
mp <- fread(file.path(RAW, "amphib_only_map.txt"))[`#SampleID` %in% colnames(R)]; setnames(mp, "#SampleID", "SampleID")
mp <- mp[match(colnames(R), SampleID)][, hsp := paste(tolower(trimws(host_genus)), tolower(trimws(host_species)))]
LV4 <- c(CLASS_LEVELS, "Not detected")
fun_class <- function(ni, nf, n) fifelse(ni / n >= 0.5, "Inhibitory", fifelse(nf / n >= 0.5, "Enhancing", "No effect / mixed"))
E <- X[, .(sOTU = unlist(strsplit(sOTUs_3mm, ";"))), by = .(FastaID, Bd_function, inhib_score, host_norm)]

## taxon-level units for the global and when-present definitions
U1 <- S[, .(unit = sOTU, definition = "Global mean", cls = abund_class(mean_rel_pct), function_class, median_score)]
U2 <- S[, .(unit = sOTU, definition = "Mean when present", cls = abund_class(mean_when_present_pct), function_class, median_score)]
## host-specific units
spN <- mp[, .N, by = hsp]; okh <- spN[N >= 5 & !grepl(" (sp|aff|cf)$", hsp), hsp]
E3 <- E[host_norm %in% okh]
hm <- sapply(okh, function(h) Matrix::rowMeans(R[, mp$hsp == h, drop = FALSE]))
U3 <- E3[, .(n_isolates = uniqueN(FastaID), n_inhib = sum(Bd_function == "Inhibitory"), n_fac = sum(Bd_function == "Facilitating"),
             median_score = median(inhib_score, na.rm = TRUE)), by = .(sOTU, host_norm)]
U3[, host_mean_pct := hm[cbind(match(sOTU, rownames(hm)), match(host_norm, colnames(hm)))]]
U3[, cls := factor(fifelse(host_mean_pct == 0, "Not detected", as.character(abund_class(host_mean_pct))), levels = LV4)]
U3[, function_class := fun_class(n_inhib, n_fac, n_isolates)]
fwrite(U3, file.path(OUT, "host_specific_units.csv"))
message("host-specific: ", uniqueN(E3$FastaID), " isolates, ", length(unique(E3$host_norm)), " host species, ", nrow(U3), " units")
U <- rbind(U1, U2, U3[, .(unit = paste(sOTU, host_norm), definition = "Host-specific (source host)", cls, function_class, median_score)])
U[, cls := factor(as.character(cls), levels = LV4)]
U[, definition := factor(definition, levels = c("Global mean","Mean when present","Host-specific (source host)"))]
summ <- U[, .(n_units = .N, pct_inhibitory = 100 * mean(function_class == "Inhibitory"), pct_no_effect = 100 * mean(function_class == "No effect / mixed"),
              pct_enhancing = 100 * mean(function_class == "Enhancing"), n_scored = sum(!is.na(median_score)), median_score = median(median_score, na.rm = TRUE)),
          by = .(definition, cls)][order(definition, cls)]
tests <- U[cls != "Not detected", .(chi2_inhibitory_p = chisq.test(table(droplevels(cls), function_class == "Inhibitory"))$p.value,
                                    kruskal_score_p = kruskal.test(median_score ~ cls)$p.value), by = definition]
fwrite(summ, file.path(OUT, "abundance_class_summary.csv")); fwrite(tests, file.path(OUT, "abundance_class_tests.csv"))
print(summ, digits = 3); print(tests, digits = 3)

## plotting helpers
xl <- c("Dominant (>=1%)" = "Dom.\n\u22651%", "Intermediate (0.1-1%)" = "Interm.\n0.1\u20131%", "Rare (<0.1%)" = "Rare\n<0.1%", "Not detected" = "Not\ndet.")
row_panels <- function(D, ylab_n, ttl, t) {
  n <- D[, .N, by = cls]
  a <- ggplot(n, aes(cls, N)) + geom_col(fill = "grey30", width = 0.6) + geom_text(aes(label = N), vjust = -0.4, size = 2) +
    scale_x_discrete(labels = xl, drop = TRUE) + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + labs(x = NULL, y = ylab_n, title = ttl) + theme_pub()
  f <- D[, .N, by = .(cls, function_class)][, pct := 100 * N / sum(N), by = cls]
  f[, function_class := factor(function_class, levels = rev(names(COL_FUN)))]
  b <- ggplot(f, aes(cls, pct, fill = function_class)) + geom_col(width = 0.6, colour = "white", linewidth = 0.2) +
    geom_text(data = f[pct >= 8], aes(label = round(pct)), position = position_stack(vjust = 0.5), size = 2, colour = "white") +
    scale_fill_manual(values = COL_FUN, breaks = names(COL_FUN)) + scale_x_discrete(labels = xl) +
    annotate("text", x = (nlevels(droplevels(D$cls)) + 1) / 2, y = 112, size = 2, label = sprintf("inhibitory vs not: \u03c7\u00b2 p = %.2g", t$chi2_inhibitory_p)) +
    coord_cartesian(ylim = c(0, 118)) + labs(x = NULL, y = "Units by Bd function (%)") + theme_pub() + theme(legend.position = "top")
  c <- ggplot(D[!is.na(median_score)], aes(cls, pmax(pmin(median_score, 145), -195))) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
    geom_jitter(width = 0.17, size = 0.4, alpha = 0.35, colour = COL_FUN["Inhibitory"]) + geom_boxplot(outlier.shape = NA, fill = NA, width = 0.55) +
    annotate("text", x = Inf, y = 175, hjust = 1.05, size = 2, label = sprintf("KW p = %.2g", t$kruskal_score_p)) +
    scale_x_discrete(labels = xl) + coord_cartesian(ylim = c(-200, 180)) + labs(x = NULL, y = "Inhibition score\n(median, %)") + theme_pub()
  list(a, b, c)
}
defs <- levels(U$definition)
pl <- unlist(lapply(defs, function(d) row_panels(U[definition == d], if (d == defs[3]) "Taxon \u00d7 host\nunits (n)" else "Bd-tested taxa (n)", d, tests[definition == d])), recursive = FALSE)
leg <- ggplot(data.table(f = factor(names(COL_FUN), levels = names(COL_FUN)), y = 1), aes(f, y, fill = f)) + geom_col() + scale_fill_manual(values = COL_FUN) + theme_void() + theme(legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 7))
p <- wrap_plots(pl, ncol = 3, widths = c(0.85, 1, 1.1)) + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect") & theme(legend.position = "top")
ggsave(file.path(FIG, "abundance_class_definitions.png"), p, width = 7.4, height = 7.6, dpi = 300)
p1 <- wrap_plots(row_panels(U[definition == "Global mean"], "Bd-tested taxa (n)", "Taxa (global classes)", tests[definition == "Global mean"]), ncol = 3, widths = c(0.85, 1, 1.1)) + plot_annotation(tag_levels = "a")
ggsave(file.path(FIG, "abundance_class_function_taxa.png"), p1, width = 7.4, height = 2.9, dpi = 300)

## isolate-level version (global classes; isolate inherits class of its sOTU)
I <- merge(E, S[, .(sOTU, cls = abund_class(mean_rel_pct))], by = "sOTU")[, .SD[1], by = FastaID]
I[, function_class := fifelse(Bd_function == "Inhibitory", "Inhibitory", fifelse(Bd_function == "Facilitating", "Enhancing", "No effect / mixed"))]
I[, median_score := inhib_score]
ti <- data.table(chi2_inhibitory_p = chisq.test(table(I$cls, I$function_class == "Inhibitory"))$p.value, kruskal_score_p = kruskal.test(median_score ~ cls, data = I)$p.value)
pi <- row_panels(I, "Bd-tested isolates (n)", "Isolates (global classes)", ti); pi[[3]] <- pi[[3]] + labs(y = "Inhibition score\n(per isolate, %)")
ggsave(file.path(FIG, "abundance_class_function_isolates.png"), wrap_plots(pi, ncol = 3, widths = c(0.85, 1, 1.1)) + plot_annotation(tag_levels = "a"), width = 7.4, height = 2.9, dpi = 300)
fwrite(I[, .(n_isolates = .N, pct_inhibitory = 100 * mean(function_class == "Inhibitory"), pct_enhancing = 100 * mean(function_class == "Enhancing"),
             n_scored = sum(!is.na(median_score)), median_score = median(median_score, na.rm = TRUE)), by = cls][order(cls)], file.path(OUT, "abundance_class_summary_isolates.csv"))
print(fread(file.path(OUT, "abundance_class_summary_isolates.csv")), digits = 3)
message("NOTE: isolate-level p-values are inflated by non-independence (isolates sharing an sOTU/study).")
