# Bd-inhibitory bullfrog skin isolates: analysis bundle (2026-11-02)

Start with **project_summary.md** (repository root). It summarises every analysis and its results, with relative links to the files in this bundle.

## Isolate IDs
All files use the IDs BF1–BF7. **isolate_ID_crosswalk.csv** maps them to the original labels (78E, BF170, BF261, the CSWG105 cryobin, CSWG105, M1241, sfa6.6).
- **Where the renaming was applied:** table contents, workbook sheets, scripts, the summary, file names, and the file names and contents inside the antiSMASH result zips.
- **Redrawn figures:** the two that carried isolate labels, antiSMASH_BGC_stacked_bar.png and AmphiBac_rank_abundance.png.
- **Old genus nicknames removed:** Starkeya, Kineo, Renibacterium and Microterricola no longer appear as isolate labels. These names still appear where they are real genus names in reference or BLAST data (e.g. *Starkeya novella*, *Renibacterium salmoninarum*).
- **First MiGA batch dropped:** the 9-report batch (HB58, BF246, and the earlier *Bacillus velezensis* isolate also labelled BF4) is not in this export, to avoid two isolates sharing the ID BF4.
- **Session log:** this is a historical record of the code that ran. Its cells that parsed that batch still name HB58 and BF246, and they call the earlier *Bacillus* isolate "BF4old".
- **AmphiBac data:** the public AmphiBac files and the AmphiBac-derived match tables were not altered. Some AmphiBac isolates of other studies have names containing these strings (e.g. Panama78E).
- **Project store:** artifact_index.csv lists files under their new names, but the stored project copies keep their original names.

## Contents
| Folder | What it holds |
|---|---|
| project_summary.md | Methods, results, decisions and open items for every analysis (sections 1–13) |
| outputs/01_MiGA_taxonomy | MiGA TypeMat summary of the 7 genomes |
| outputs/02_antiSMASH | BGC summary workbook, region table, stacked-bar figure, and the full antiSMASH 7.1.0 result zip per genome (open index.html) |
| outputs/03_genus_abundance | Target-genus abundance: meta-analysis, bullfrog, salamander (summary workbooks and per-sample tables) |
| outputs/04_AmphiBac_vs_meta | AmphiBac 2025.1-tested vs meta-analysis: matching, cultured share, rank abundance, source host, abundance class × Bd function (3 definitions); isolate→sOTU candidate table (all tied sOTUs kept) and per-sOTU abundance stats |
| outputs/05_multihit_sensitivity | Sensitivity of the abundance-class results to how isolates tied between sOTUs and sOTUs with conflicting labels are resolved |
| outputs/06_within_sOTU_Bd_variance | Within-sOTU spread of Bd inhibition (current, post-label-filter versions) |
| outputs/07_panelA_separate_manuscript | Stand-alone panel A package (data, R script, methods) for the other manuscript |
| outputs/08_functional_genomics | CAZyme/mucin, KEGG module and comparative-genomics results |
| outputs/09_hormone_screen | Steroid and hormone-related KO tables |
| outputs/10_conditional_rarity | Consistently vs conditionally rare sOTUs and their taxonomy |
| code/amphib_R | R pipeline (scripts 00–08, run_all.R, BLAST script); see its README for the input files to add |
| code/multihit_sensitivity.py | Stand-alone Python script for outputs/05; verified to reproduce them from the AmphiBac files and the meta-analysis OTU table |
| code/sOTU_Bd_panelA | The panel-A package, unzipped |
| code/session_log | Every successful code cell from the analysis session (181 cells), numbered in run order, with index.csv giving the time, environment and files written. This is the full record of the Python, bash and R work, including genome, antiSMASH, MiGA and salamander steps that are not in the R pipeline. `<ARTIFACT name>` marks a project file that was read; `host.*` calls are platform functions for reading and saving files. The 13 cells that errored are not included. |
| code/functional_genomics | Scripts for outputs/08 |
| code/hormone_KO_screen | Script for outputs/09 (run end to end from the KofamScan detail output) |
| code/conditional_rarity | Scripts for outputs/10 (run end to end from the meta-analysis files) |
| artifact_index.csv | Every project file with its artifact ID, its role (output / your input / public input) and its location in this bundle |

## Not included
Your own input data, listed in artifact_index.csv as "your input (not bundled)":
- the meta-analysis OTU table, map and taxonomy;
- the bullfrog and salamander ASV tables and FASTAs;
- the 7 genome assemblies;
- the MiGA PDF reports (BF1–BF7 under their original file names; HB58, BF4 *Bacillus* and BF246 from the dropped batch).

The R pipeline and the multi-hit script expect the meta-analysis and bullfrog files in their input folders, as described in their READMEs.

## Software
- Python 3.11: pandas, numpy, scipy, matplotlib, openpyxl.
- R 4.5: data.table, Matrix, ggplot2, patchwork, scales.
- BLAST+ 2.17.0, against the NCBI 16S_ribosomal_RNA database.
- antiSMASH 7.1.0 (public web server).
- MiGA TypeMat (web).
