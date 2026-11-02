# Panel A: within-sOTU spread of Bd inhibition (stand-alone)

Run with `Rscript panelA_sOTU_Bd_spread.R` from this folder (requires R ≥4.1, data.table and ggplot2; cairo is used for the PDF).

- `data/panelA_isolate_scores.csv`: one row per scored isolate (4,551). Columns: sOTU_id (sOTU_001–174, ranked by median score), sOTU_V4_sequence (90-bp meta-analysis sOTU), FastaID (AmphiBac ID), Reference (study DOI), host_species, Bd_function, inhib_score = (1 − proportional growth) × 100, mean_rel_pct and prevalence_pct (sOTU in the meta-analysis, 2,411 samples).
- `panelA_sOTU_Bd_spread.R`: computes the per-sOTU summary and draws the panel.
- `output/`: panelA_sOTU_summary.csv, panelA_sOTU_Bd_spread.png (300 dpi) and .pdf (vector).
- `methods_panelA.md`: methods paragraph, numbers to cite and a suggested legend.

The input table was derived from AmphiBac 2025.1-tested and the meta-analysis OTU table by scripts 01, 04 and 05 of the companion `amphib_R` project. The isolate-to-sOTU assignment and the filters are described in the methods.
