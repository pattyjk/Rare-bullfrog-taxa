# amphib_R — R pipeline for the meta-analysis, bullfrog and AmphiBac analyses

Reproduces the target-genus abundance (meta-analysis and bullfrog) and all AmphiBac-vs-meta-analysis analyses in R. Genome (MiGA, antiSMASH) and salamander analyses are not included, by design.

## Setup
1. Install R ≥ 4.3 with the packages data.table, Matrix, ggplot2, patchwork and scales.
2. Put your own input files into `data/raw/` under these names (they are not bundled):
   - `amphib_only_otu_table.zip` (meta-analysis sOTU × sample table)
   - `amphib_only_map.txt` (meta-analysis sample map)
   - `meta_analysis_taxonomy.txt` (meta-analysis taxonomy)
   - `bullfrog_asv_table.txt` (bullfrog ASV table)
   - `bullfrog_dna-sequences.fasta` (bullfrog 120-bp ASV sequences)
3. Already provided in `data/raw/`:
   - AmphiBac 2025.1-tested FASTA and metadata (public, from github.com/AmphiBac/AmphiBac-Database);
   - `ncbi_refs/*.fasta`, the NCBI GenBank 16S reference sequences for 11 genera (downloaded 2026-10-16);
   - `bullfrog_asv_blast_16SRefSeq.tsv`, the BLAST output, re-created on 2026-11-02 with `blast_bullfrog_asvs.sh`.

## Run
From this folder: `Rscript run_all.R`. It runs the scripts in order and writes `output/sessionInfo.txt`. Peak memory is about 5 GB.
- Scripts 01 and 04 write `data/derived/`, which the later scripts read.
- `rel_abund_sparse.rds` (17 MB) is not bundled; script 01 regenerates it.

| Script | What it does | Main outputs |
|---|---|---|
| 00_functions.R | Paths, thresholds, helper functions (V4 window after 515F, abundance class, Bd-label harmonisation) | — |
| 01_prepare_meta.R | Parses the OTU table, keeps samples with ≥1,000 reads (2,411), computes per-sOTU mean, prevalence, mean-when-present and rank | data/derived/sOTU_stats.csv, sample_totals.csv |
| 02_meta_genus_abundance.R | Target-genus prevalence and abundance in the meta-analysis (meta taxonomy at 97%, exact V4 match to NCBI references, ≥2-reference and genus-unique variants) | output/meta_genus_* |
| 03_bullfrog_genus_abundance.R | The same for the bullfrog data, plus BLAST top and tied hits | output/bullfrog_* |
| 04_amphibac_match.R | Matches AmphiBac isolates to sOTUs: exact containment, else nearest ≤3 mismatches over the V4 window; ties go to the most abundant sOTU | data/derived/amphibac_isolate_matches.csv |
| 05_sOTU_function_variance.R | Within-sOTU spread of Bd inhibition scores; ICC variance partitioning | output/sOTU_Bd_function_variance.csv, Bd_variance_partitioning.csv, figures/sOTU_Bd_variance.png |
| 06_cultured_share.R | Share of sOTUs and reads with a cultured match, by abundance class and per sample | output/cultured_by_abundance_class.csv, per_sample_cultured_share.csv, figures/rank_abundance.png, cultured_vs_abundance.png |
| 07_abundance_class_function.R | Bd function by abundance class under three definitions (global, mean-when-present, host-specific) | output/abundance_class_*, figures/abundance_class_*.png |
| 08_source_host.R | Abundance on the isolate's source host vs other hosts | output/source_host_*, figures/source_host_abundance.png |
| optional_fetch_ncbi_refs.R | Optional: re-downloads the NCBI references (untested here; E-utilities timed out) | — |

## Key settings
- Abundance classes use mean relative abundance across all 2,411 samples: dominant ≥1%, intermediate 0.1–1%, rare <0.1%.
- Inhibition score = (1 − Proportional.Growth.Bd) × 100. It exists only for 96-well cell-free supernatant (CFS) assays.
- AmphiBac matching uses ≤3 mismatches over 90 bp (allowing for Sanger error). Ambiguity codes count as mismatches. Only sOTUs with reads are used.
- Tied nearest sOTUs go to the most abundant one. This choice matters for the abundance-class χ² test; see ../multihit_sensitivity.py and outputs/05_multihit_sensitivity.

## Provenance and checks
This folder was rebuilt on 2026-11-02 by replaying the recorded session code (see ../session_log, cells 136–153). Its outputs match the earlier deliverables, with one exception. In `output/bullfrog_genus_summary.csv`, *Arthrobacter* under `exact120_unique_min2refs` has 0 ASVs here vs 1 in the earlier Python workbook. The ASV's exact match is shared with *Paeniglutamicibacter*, so 0 is consistent with the "genus-unique" definition.
The NCBI 16S BLAST database used for the re-run (29 Oct 2026, 27,646 sequences) is one sequence smaller than the original (27,647); the BLAST-based results are identical.
