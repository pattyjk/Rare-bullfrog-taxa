# KEGG Orthology annotation and module completeness (BF1–BF7)

## Inputs
Prokka protein sets (`fg_inputs/<ID>/<ID>.faa`). Counts: BF1 4,677; BF2 3,775; BF3 5,608; BF4 3,888; BF5 3,993; BF6 3,601; BF7 5,567 proteins (31,109 total). All proteins were used. Protein IDs were prefixed `<ID>__` and the seven sets were run together in one KofamScan job. KOfam thresholds are per profile, so this does not change any assignment.

## Tools and databases
- KofamScan `exec_annotation` 1.3.0 (Ruby 3.4.2), HMMER 3.4 (Aug 2023), GNU parallel 20260922, Python 3.11.17 / pandas.
- KOfam (https://www.genome.jp/ftp/db/kofam/): `profiles.tar.gz` (HTTP Last-Modified 27 Oct 2026; 27,981 KO profiles) and `ko_list.gz` (Last-Modified 26 Oct 2026; 28,504 KOs). The full profile set was used, not `prokaryote.hal`. The database is not redistributed.
- KEGG modules: KEGG REST `list/module` and `get/<M>`, retrieved 2026-11-02 (`info/module` reports 583 entries, "Last update 2026/10/28"). Snapshot: `kegg_module_definitions_snapshot_2026-11-02.tar.gz`. It contains the raw entries plus the KO lists for pathway maps map00984, map00362, map00627, map00622, map00623, map00624, map00626, map00621, map00633, map00643, map00364, map00361, map00791, map00930 and map02020 (KEGG REST `link/ko/<map>`, same date).
- CPU: `--cpu 3` gives GNU parallel `-j 3` with `hmmsearch --cpu 1` per profile. Wall time was 40 min; user CPU 115 min, about 2.9 cores on average.

## Thresholds
- **KO assignment (primary):** only significant KofamScan hits (`*`), i.e. HMMER score ≥ the KOfam adaptive per-KO threshold (full or domain score, as set in ko_list). All significant KOs per protein are kept (`n_sig_KOs_for_protein`). KO counts = number of distinct proteins assigned that KO.
- **Module completeness:** module DEFINITIONs were parsed with `module_completeness.py`, which has its own parser. Precedence: `+`/`-` complex components > space (AND, sequential steps) > `,` (OR, alternatives). Parentheses group terms.
  - `-K` marks an optional component and is ignored.
  - Steps that are entirely optional, and `--` steps (no known KO), are excluded from the denominator.
  - Embedded module IDs count as satisfied only if that module is complete.
  - Multi-line definitions are joined as additional steps.
  - A step is present if any alternative has all of its non-optional KOs.
  - completeness = present steps / scorable top-level steps. If the top level is itself an OR of paths, the best path is used.
  - `complete` = 1.0. `ge075` = ≥0.75.
- **Curated features:** a feature with alternative modules (e.g. lysine pathways) takes the maximum completeness among them. Possible amino-acid auxotrophy = best completeness < 0.5.
- **Urease structural call:** ureC (K01428) plus [ureB K01429 and ureA K01430, or fused ureAB K14048].
- **Two-component systems:** a pair counts as present when both the sensor-kinase KO and the cognate response-regulator KO are significant.
- **Sensitivity only (not a primary call):** `sensitivity_completeness_relaxed` adds non-significant hits with score ≥ 0.5 × the adaptive threshold and E ≤ 1e-10. This adds 15,456 isolate–KO pairs, many of them paralog cross-hits, so it is an upper bound. `curated_missing_steps_subthreshold_hits.tsv` lists, for every missing step of a curated module, the best sub-threshold hit per KO (score, threshold, ratio, E-value).

## Caveats
- KO presence shows genetic capacity only, not expression or activity. Absence at the KOfam threshold is not proof of absence. Divergent homologs (common in novel genera/species such as BF2–BF6) can fall below adaptive thresholds. Examples:
  - BF4 ureB LJHHHOOL_00331 scores 147.6 against a 172.7 threshold, but Prokka annotates a complete ureDGFCBA operon.
  - BF1 prephenate dehydratase BEDMLJNG_02084 scores 324.6 against 326.2.
  - The monomeric isocitrate dehydrogenase of BF5/BF7 (K00031) scores 106/94 against 153.6.
  - BF3/BF7 cysteine synthase K01738 scores 397/400 against 423.
- KEGG module definitions do not cover every lineage's variant enzymes. The Phe/Tyr final aminotransferase step (K00832/K00838) is unassigned in all isolates. Broad-specificity aminotransferases (K00812, K00817, K14267) are present, so <0.5 Phe/Tyr scores do not by themselves indicate auxotrophy. Similar issues: ThiS/ThiF/ThiI/ThiH sulfur-transfer steps of thiamine synthesis; actinobacterial APS reductase in assimilatory sulfate reduction (M00176 expects APS kinase K00860); the ubiquinone vs menaquinone split (Alphaproteobacteria BF1/BF2 vs Actinomycetota).
- BF4 and BF5 are two bins from one mixed culture. Binning errors could move genes between them.
- Draft assemblies: genes split across contig ends may be missed.
