# session cell 170 | 2026-10-28 18:36 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/project_summary.md","sha256":"c5c5ef326d0476ba66ab9d2c14aacd926f26f315fb24f80bc5c6d592bdfd113a"}]

A=host.artifact_marker
L=lambda name,vid: f"[{name}]({A(vid)})"
txt=f"""# Project summary: Bd-inhibitory bullfrog skin isolates as rare / novel taxa
*Status as of 2026-10-18. All numbers come from the saved artifacts linked below.*

## 1. Goal
Bacteria cultured from bullfrog (*Lithobates catesbeianus*) skin inhibit *Batrachochytrium dendrobatidis* (Bd) at about 100% in vitro. The aim is to show that these isolates (i) are undescribed or poorly known taxa, (ii) are rare or absent in amphibian skin amplicon surveys, and (iii) carry novel biosynthetic potential. They are set against the background of the AmphiBac isolate database and the amphibian skin microbiome meta-analysis.

## 2. Genome taxonomy (MiGA, TypeMat)
Tables: {L('MiGA_isolate_summary_7.xlsx','3c39667b-c730-4221-b772-eecaee89a9b4')} (7 genomes) and {L('MiGA_isolate_summary.xlsx','4914d5c2-451b-42ce-ada1-f2570b2b22f0')} (first 9 reports). All genomes were rated Excellent quality by MiGA (100% complete, ≤4.3% contamination).

| Isolate | Closest type genome | Metric | Interpretation |
|---|---|---|---|
| BF7 | *Rhodococcoides fascians* DSM 43673 | ANI 98.05% | *R. fascians* (≥95% ANI; MiGA's p-value call says "novel species") |
| BF3 | *Kineococcus rhizosphaerae* DSM 19711 | AAI 74.21% | *Kineococcus* sp. nov. |
| BF6 | *Psychromicrobium lacuslunae* IHBB 11108 | ANI 83.30% | *Psychromicrobium* sp. nov. |
| BF5 | *Arthrobacter humicola* JCM 15921 | ANI 88.16% | *Arthrobacter* sp. nov. |
| BF4 | *Cryobacterium soli* GCJ02 | ANI 90.97% | *Cryobacterium* sp. nov.; second organism in the BF5 culture (different family); no 16S gene recovered |
| BF2 | *Ancylobacter pratisalsi* DSM 102029 | AAI 70.15% | Xanthobacteraceae, possibly a new genus |
| BF1 | *Sphingomonas paucimobilis* | ANI 99.49% | *S. paucimobilis*; contaminant *Bacillus* 16S on one low-coverage contig (NODE_103) |

Isolates covered only in the first batch: HB58 is *Acinetobacter oleivorans* (ANI 95.25%), BF4old is *Bacillus velezensis*, and BF246 is *Microbacterium* sp. nov.
Taxonomy note: *Rhodococcoides* is the recent split of *R. fascians*, *R. corynebacterioides*, *R. kroppenstedtii*, *R. kyotonensis*, *R. navarretei*, *R. trifolii*, *R. yunnanensis* and *R. luteus* from *Rhodococcus*. Older 16S databases still label these as *Rhodococcus*.

## 3. Biosynthetic potential (antiSMASH 7.1.0, relaxed)
Files: {L('antiSMASH_BGC_summary.xlsx','874db76f-e23b-47bd-8a49-017504f23ea1')}, {L('antiSMASH_regions.csv','107585dd-e41f-4250-9007-7b68bb9b313f')}, and one full-result zip per genome.
- There are 57 regions in total: BF7 23, BF3 9, BF6 9, BF5 6, BF4 4, BF2 4, BF1 2.
- Only 3 regions reach ≥75% similarity to a known MIBiG cluster: ε-poly-L-lysine in BF7 and in the cryobin, and zeaxanthin in BF1. Most regions are below 50%, i.e. potentially novel.
- The two BF5 assemblies share NAPAA, T3PKS and betalactone region types.

![BGC types per genome]({A('3216d9f7-d15c-401c-8371-7242e5a4bd11')})

## 4. Target genera in amplicon data
### 4a. Meta-analysis
2,411 samples with ≥1,000 reads. Methods: BLAST (≥97%) and exact V4 match against NCBI 16S references for 11 genera. Files: {L('amphibian_genus_abundance.xlsx','487b5319-6b32-4e83-9134-dc48b16f5f10')}, per-sample {L('amphibian_genus_relabund_by_sample.csv','7ce1bc36-0aaa-40ec-b461-d22e24b98de1')}.
- *Sphingomonas*: present in 80.9% of samples (mean 0.78%). *Arthrobacter* 41.7%, *Rhodococcus* 32.4%, *Kineococcus* 14.4%, *Ancylobacter* 13.7%.
- *Rhodococcoides* and *Psychromicrobium*: no sOTU has a BLAST top hit to either genus. The exact V4 matches that exist are shared with *Rhodococcus* or *Arthrobacter*.
- *Cryobacterium*: exact-V4 prevalence is 23.9%, falling to 16.9% when at least 2 supporting reference records are required.

### 4b. Your bullfrog dataset
42 samples, 3,700 ASVs of 120 bp; P10F3 and P4F1 were dropped for having <1,000 reads. Files: {L('bullfrog_genus_abundance.xlsx','9836893a-0a4b-41bb-ba47-0c6ae73c7590')}, {L('ASV_target_genus_matches.xlsx','a2be2d04-8839-402f-997d-b9252c2764f3')}.
- *Sphingomonas*: 0.93–1.42%, in 100% of samples. *Ancylobacter* about 0.08%.
- *Rhodococcoides* candidate (ASV 5cb28bd0; tied with *Antrihabitans*/*Rhodococcus*): 0.04%, in 50% of samples.
- *Kineococcus* 0.01–0.03%; *Arthrobacter* 0.008%.
- Not detected: *Psychromicrobium*. The *Cryobacterium* signal rests only on single-record matches and is likely spurious.

### 4c. Red-backed salamander (RBS1–15)
File: {L('salamander_genus_abundance.xlsx','dd16e100-c8c7-4cb2-b88e-cf5ea8204584')}.
- *Sphingomonas* 0.80% (part of this signal has BLAST top hits to *Parasphingomonas*). *Ancylobacter* 0.11%, in 53% of samples.
- Not detected: *Rhodococcoides*, *Psychromicrobium* and *Kineococcus*.

### 4d. Your isolates' own V4 sequences (from the genomes) in the meta-analysis
- BF5: 0.22% (rank 55). BF1: 0.052%. BF3: 0.012%. BF6: 0.0007%.
- Absent: BF7 and BF2.
- AmphiBac isolates that share BF7's V4 are mostly non-inhibitory (19 of 28). All 3 that share BF3's V4 are non-inhibitory.

## 5. AmphiBac 2025.1 (tested) vs the meta-analysis
Isolates were matched to sOTUs by exact containment of the 90-bp sOTU in the isolate 16S, otherwise by ≤3 mismatches over the V4 window (allowing for Sanger error). Files: {L('AmphiBac_meta_abundance.xlsx','33c8c7ee-0678-4428-92b2-fc3490b532e1')}, {L('AmphiBac_near_matches.xlsx','02bd023d-f792-4d52-84d4-34c94e8f3920')}.
- **Matching:** 6,610 of the 6,729 V4-spanning isolates match an sOTU with reads. 1,144 isolates lack V4 coverage (median length 488 bp).
- **Cultured share by abundance class:** 15/15 dominant, 55/82 intermediate, 494/48.3k rare sOTUs have a cultured representative. Per sample, cultured sOTUs are a median 20% of sOTUs but 68% of reads.
- **Detected vs undetected isolates:** inhibition rates are similar, 49% vs 43% (p = 0.32).

![Rank abundance]({A('9b827439-f9ce-496c-a860-bf3339be6a37')})
![Cultured vs abundance]({A('220050de-8253-43e7-8e00-2c1e2a620e6f')})

**Source host** ({L('AmphiBac_source_host_abundance.xlsx','1ce89640-3567-47af-87cb-80a6f015cea3')}):
- Isolates are usually less abundant on the host species they were cultured from than elsewhere: the source host is higher in only 26% of the 1,626 sOTU × host units (pooled).
- Host by host (1,476 units), the source host ranks above a median 54% of other host species.
- The taxon is absent on its source host in 35% of units, and the source host is the top host only 40 times.

![Source host]({A('5e4f0031-2b62-4bca-861f-e61d7cbb5c07')})

## 6. Abundance class vs Bd function
Classes: dominant ≥1%, intermediate 0.1–1%, rare <0.1%. Files: {L('AmphiBac_abundance_class_function.xlsx','9fe895e2-f783-47fb-a6fe-01cb90cbfb07')}, {L('AmphiBac_abundance_class_definitions.xlsx','c6b6f3c7-566c-449f-bcb0-c9d33157accb')}.

| Definition | n (D / I / R) | % inhibitory (D / I / R) | p (χ²) |
|---|---|---|---|
| Global mean, taxa | 15 / 55 / 472 | 67 / 49 / 39 | 0.048 |
| Global mean, isolates | 1,017 / 1,797 / 3,158 | 68 / 52 / 40 | (inflated, non-independent) |
| Mean when present, taxa | 65 / 278 / 199 | 57 / 44 / 33 | 0.0013 |
| Host-specific (59 hosts; + 523 not detected) | 88 / 230 / 638 | 62 / 55 / 44 (not detected 42) | 0.0004 |

**Conclusion:** under every definition, rare taxa are *less* often Bd-inhibitory, not more. Framing for the paper: your isolates are an exception to the database-wide trend. They are rare or undetected yet strongly inhibitory, which shows that inhibition is a strain-level trait not predictable from 16S.

![Abundance class vs function]({A('ba58bb37-446c-4c19-a534-3510f43ed1c0')})
![Alternative definitions]({A('0f3418be-68ff-42f3-b136-11d8128ed69e')})

## 7. Variation in Bd inhibition within each sOTU
Score = (1 − proportional Bd growth) × 100, from 96-well cell-free supernatant (CFS) assays only. Current values, after restricting to the three Bd labels:
- **Spread:** 174 sOTUs have ≥5 scored isolates (4,551 isolates, 16 studies). The median within-sOTU SD is 39.6 points.
- **Opposite effects:** 74.7% of sOTUs contain both an isolate scoring >50 and one scoring <0 (enhancing Bd growth).
- **No class effect:** median SD by class is 37.0 / 38.8 / 40.6 (Kruskal–Wallis p = 0.92).
- **Label agreement** (sOTUs with ≥3 isolates): a median of 67% of isolates share the majority label, and 13% of sOTUs are unanimous.
- **Variance partitioning:** sOTU explains 0.15 of the variance (ICC), study 0.23, and sOTU × study 0.39. 30% of the 146 multi-study sOTUs show a significant study effect.

The older saved files {L('sOTU_Bd_variance.png','fec1c944-eada-4a13-9160-3791c6b5b46f')}, {L('sOTU_Bd_function_variance.csv','078a8a73-8ddf-43b5-b701-0e011a472d1f')} and {L('Bd_variance_partitioning.csv','5d117ff4-ebbe-4577-a119-88b601b49707')} predate this filter. They show 175 sOTUs, 11% unanimous, and 29% of 147 sOTUs.

**Panel A, packaged for the separate manuscript:** {L('sOTU_Bd_panelA.zip','c9243529-ebfa-4c5f-941f-c4ef3fb1cdb9')} (data, R script, figure, methods; abundance colouring removed).
![Panel A]({A('928aeecd-23f0-4a4b-a1d7-ee91c5fba0f9')})

## 8. Decisions recorded
- AmphiBac matching uses ≤3 mismatches over V4 (your choice, to allow for Sanger error).
- Taxon-level p-values are the conservative ones; isolate-level tests are inflated by non-independence.
- Genome and salamander analyses are kept out of the R package.

## 9. Open items
1. **R project (`amphib_R`)**: scripts 01–08 had been written and checked against the Python results (script 08 reproduced 1,626 units, 54%, 35% and 40). However, the project was never saved as an artifact and was lost when the workspace was cleaned. It needs to be rebuilt and saved as a zip with run_all.R, a README and the BLAST script.
2. The optional NCBI re-download script timed out against NCBI E-utilities and is untested. The reference sequences used are saved in {L('ref16S_11genera_ncbi.tar.gz','05d759da-dbb2-4d8a-84b1-fb5fc6866353')}.
3. Offered but not run: the ≥2-record *Cryobacterium* filter for the per-sample meta-analysis files, and bullfrog-only abundance classes.
"""
open('project_summary.md','w').write(txt); print(len(txt))