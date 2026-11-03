# CAZyme, mucin O-glycan and fungal cell-wall screen (BF1–BF7): methods

## Inputs
Prokka proteomes (`fg_inputs/<ID>/<ID>.faa`; 31,109 proteins in total: BF1 4,677; BF2 3,775; BF3 5,608; BF4 3,888; BF5 3,993; BF6 3,601; BF7 5,567) and Prokka GFF3 files for the CDS coordinates. No proteins were removed (BF1 contaminant contig NODE_103 carries no CDS). BF4 and BF5 are two bins from one mixed culture, so a gene could be assigned to the wrong bin.

## Tools and databases
- pyhmmer 0.12.3 (HMMER3 API), Python 3.12, pandas 3.0.6 (conda env `chytrid-cazy`). All searches used `cpus=2`.
- dbCAN HMMdb V13 (826 HMMs; https://pro.unl.edu/dbCAN2/download/Databases/V13/dbCAN-HMMdb-V13.txt), downloaded 2026-11-02.
- Pfam and NCBIfam/TIGRFAM HMMs (67 models), downloaded from the InterPro API (`/entry/{pfam|ncbifam}/<acc>?annotation=hmm`) on 2026-11-02. Accession versions are listed in `panel_hmm_versions.tsv`.
- DeepSig (Savojardo et al. 2018, *Bioinformatics* 34:1690), GitHub BolognaBiocomp/deepsig commit 727941b, run under TensorFlow/Keras 2.15.1 (conda env `deepsig`, Python 3.10).

## Thresholds
- **dbCAN** (same as the earlier run_cazome.py/run_hmm.py): hmmsearch E ≤ 1e-5, domE ≤ 1e-5. A domain is kept if its i-Evalue is < 1e-15 when the HMM alignment is longer than 80 aa, otherwise < 1e-5, and its HMM coverage is > 0.35.
- **Overlap resolution** (full CAZome): within a protein, when two domain envelopes overlap by more than 50% of the shorter one, only the hit with the best i-Evalue is kept. Subfamily models (GH5_x, GH13_x, GH16_x, GT2_* …) are collapsed to the family (`substrate_map.family`).
- **Counting**: counts are the number of distinct proteins per isolate with at least one kept domain of the family. "Degradative CAZyme proteins" are distinct proteins with any GH, PL, CE, AA or CBM domain. Values are normalised per 1,000 Prokka proteins.
- **Pfam** models used gathering (GA) cutoffs. **NCBIfam/TIGRFAM** models used trusted (TC) cutoffs.
- **Operon neighbourhood**: core genes count as co-located when they are on the same contig and within ±10,000 bp of a core-gene anchor (measured from the anchor CDS edges).
- **Signal peptide**: DeepSig DNN class "SignalPeptide", using DeepSig's own thresholds (gramn 0.7575 for BF1/BF2; gramp 0.72 for BF3–BF7). Tat motif: a hit to TIGR01409 or PF10518 (TC/GA) that starts within residues 1–60.
- **Most-variable families**: ranked by log2((max+1)/(min+1)) of per-isolate protein counts, with ties broken by range (max − min).

## (1) Full CAZome
All 826 dbCAN V13 HMMs were searched against all 7 proteomes (`run_cazome_bf.py`, then `analyze_cazome_bf.py`). Families were assigned to a primary substrate with the earlier curated `substrate_map.py`, unchanged. Families missing from the map are labelled "Mixed / other"; this includes the recent families GH179 and GH188. GTs are reported separately as biosynthetic families and are excluded from the substrate matrix.

## (2) Mucin O-glycan panel
- **Families**: the 30-family dbCAN subset from run_hmm.py plus Pfam PF00884, PF03781, PF13402, PF10462, PF02128, PF00026, PF03629 and PF04131 (GA).
- **Primary counts** follow run_hmm.py: each subset family is counted on its own, with no overlap removal across families. The CAZome overlap-resolved count is given alongside. GH16 is reported as its own row because most GH16 enzymes are not mucin-active.
- **Catabolic operons**: genes were identified with HMMs only, not Prokka names.
  - Sialic acid: nanA (TIGR00683, NF003164), nanK (NF003461, NF047821), nanE (NF002231, PF04131).
  - L-fucose: fucI (TIGR01089, NF008220), fucK (TIGR02628), fucA (TIGR01086, NF004527, NF005302, NF005984, NF040649). Accessory genes: fucO (TIGR02638, NF007911), fucU (NF011949), fucP (TIGR00885, NF007524).
  - GlcNAc/GalNAc: nagA (TIGR00221, NF008371) and nagB (TIGR00502, NF001682/4/5, NF009022, NF041128, NagB-II NF046059). Aga/Kba accessory genes: agaS (NF007291, TIGR02815), GalNAc/GalN PTS subunits (NF007288, NF007289, NF007359, NF008508, NF040757, NF040760, NF040761), kbaY/gatY (TIGR01858, NF009375), kbaZ/gatZ (TIGR02810, NF012002, PF08013), lacC (TIGR01231, NF010033) and lacD (TIGR01232).
  - Sialic acid transporters: no specific NanT/SiaP/SatABC HMM exists in Pfam/NCBIfam. Only generic transporter families were used, and only when they lie in the same neighbourhood (TRAP: PF03480, PF06808, PF04290; MFS: PF07690, PF00083; ABC sugar uptake: PF01547, PF13416, PF00528). These are reported as accessory candidates with low specificity.
  - Calls: "complete core cluster" = all core genes co-located; "partial" = ≥2 co-located; "isolated" = one core gene only.
  - Diagnostic: for each core model, `operon_model_best_hits_diagnostic.tsv` gives the best hit at E ≤ 1e-5 with its bit score relative to TC/GA. Hits at ≥90% of the cutoff that do not pass it are listed separately as provisional near-threshold candidates; they are not used in the calls.

## (3) Fungal cell-wall panel
- **Families** (taken from the overlap-resolved CAZome):
  - Chitin: GH18, GH19, GH23 (mostly peptidoglycan lytic transglycosylases; listed because the task requested it), AA10, CBM5, CBM12 and CBM73. GH20 (chitobiase/hexosaminidase) is noted. GH46 and GH75 chitosanases were added as extras.
  - β-glucan: GH16, GH55, GH64, GH81, GH17, GH128 and GH132.
- **GH16 subfamily**: the best-scoring GH16 model of the protein. "GH16" means the family-level model scored highest, so no subfamily could be resolved.
- **Secretion**: DeepSig and the Tat HMM were run on every protein of every proteome (`deepsig_all_proteins.tsv`).

## Deviations and caveats
- **SignalP 6 not used**: SignalP 6 needs a DTU licence download, so DeepSig was used. The bioconda DeepSig package could not be installed because it pins TensorFlow 2.2, which has no build for this macOS platform. DeepSig was therefore installed from GitHub and its DNN ensemble called directly (`run_deepsig_bf.py`). On DeepSig's bundled test set (46 eukaryotic proteins) this reproduced the reference signal-peptide calls exactly (42/42, with no extra or missing calls).
- **No cleavage sites**: DeepSig's cleavage-site CRF is a Linux-only static binary, and its pure-Python CRF fallback fails (index errors). Cleavage positions are therefore not reported (`cleavage_pos` = "-"); only SP / transmembrane / other classes with probabilities are given. One Python 3 integer-division fix was made in `deepsig/crfdecoding.py` during this attempt; it does not affect the DNN classification. DeepSig does not distinguish lipoprotein (Sec/SPII) signals, so some "Transmembrane" calls may be signal anchors or lipoproteins.
- **Gene presence ≠ activity**: a domain hit shows genetic capacity, not demonstrated enzyme activity, substrate specificity or secretion.
- **Low-specificity families**: GH109, GH179 and GH188 share the Gfo/Idh/MocA NAD-dependent oxidoreductase fold. Here most of their hits are Prokka-annotated inositol dehydrogenases or similar oxidoreductases, so these families inflate GH counts and are flagged in the variability table. GH2/GH35/GH42 and GH29/GH95 are polyspecific, and GH23 is mainly peptidoglycan-active.
- **Draft assemblies**: assemblies are fragmented (127–2,456 contigs), so operon neighbourhoods can be split across contigs. Absence of co-location is therefore weaker evidence than presence.
