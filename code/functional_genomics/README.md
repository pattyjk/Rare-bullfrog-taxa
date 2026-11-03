# Functional-genomics scripts
- **Inputs:** `fg_inputs/<ID>/<ID>.faa|.gff|.fna`, i.e. the Prokka outputs in `Prokka/` (BF1–BF7). Proteins on the BF1 contaminant contig would be removed, but that contig carries none.
- **cazymes/:** dbCAN V13 + Pfam/NCBIfam search (`run_cazome_bf.py`), post-processing into matrices, mucin panel, operons and fungal-wall panel (`analyze_cazome_bf.py`), and the substrate map. `run_deepsig_bf.py` holds the DeepSig signal-peptide calls.
- **kegg/:** KofamScan run and database download (`run_kofamscan.sh`), the KEGG module-definition parser (`module_completeness.py`) and the curated tables (`kegg_postprocess.py`).
- **comparative/:** NCBI genome fetch/QC/ANI/dereplication (00), gene calling (01), fas/att/hyp screen (02), Panaroo + MMseqs2 pangenomes and isolate-specific genes (03), core trees (04), Pfam domains and islands (05), island sharing (06), and tree placement (07).
- **Figures:** the heatmaps in outputs/08_functional_genomics were drawn from the saved tables; the plotting code is in the main session log.
See ../../outputs/08_functional_genomics/methods/ for tools, versions and thresholds.
