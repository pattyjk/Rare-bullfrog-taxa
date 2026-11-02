# Project summary: Bd-inhibitory bullfrog skin isolates as rare / novel taxa
*Status as of 2026-11-02. All numbers come from the saved artifacts linked below.*

*Isolate IDs follow the BF1–BF7 scheme (crosswalk to the original labels: isolate_ID_crosswalk.csv).*

## 1. Goal
Bacteria cultured from bullfrog (*Lithobates catesbeianus*) skin inhibit *Batrachochytrium dendrobatidis* (Bd) at about 100% in vitro. The aim is to show that these isolates (i) are undescribed or poorly known taxa, (ii) are rare or absent in amphibian skin amplicon surveys, and (iii) carry novel biosynthetic potential. They are set against the background of the AmphiBac isolate database and the amphibian skin microbiome meta-analysis.

## 2. Genome taxonomy (MiGA, TypeMat)
Tables: [MiGA_isolate_summary_7.xlsx]({{artifact:3c39667b-c730-4221-b772-eecaee89a9b4}}) (7 genomes). All genomes were rated Excellent quality by MiGA (100% complete, ≤4.3% contamination).

| Isolate | Closest type genome | Metric | Interpretation |
|---|---|---|---|
| BF7 | *Rhodococcoides fascians* DSM 43673 | ANI 98.05% | *R. fascians* (≥95% ANI; MiGA's p-value call says "novel species") |
| BF3 | *Kineococcus rhizosphaerae* DSM 19711 | AAI 74.21% | *Kineococcus* sp. nov. |
| BF6 | *Psychromicrobium lacuslunae* IHBB 11108 | ANI 83.30% | *Psychromicrobium* sp. nov. |
| BF5 | *Arthrobacter humicola* JCM 15921 | ANI 88.16% | *Arthrobacter* sp. nov. |
| BF4 | *Cryobacterium soli* GCJ02 | ANI 90.97% | *Cryobacterium* sp. nov.; binned from the same culture as BF5 (different family); no 16S gene recovered |
| BF2 | *Ancylobacter pratisalsi* DSM 102029 | AAI 70.15% | Xanthobacteraceae, possibly a new genus |
| BF1 | *Sphingomonas paucimobilis* | ANI 99.49% | *S. paucimobilis*; contaminant *Bacillus* 16S on one low-coverage contig (NODE_103) |

Taxonomy note: *Rhodococcoides* is the recent split of *R. fascians*, *R. corynebacterioides*, *R. kroppenstedtii*, *R. kyotonensis*, *R. navarretei*, *R. trifolii*, *R. yunnanensis* and *R. luteus* from *Rhodococcus*. Older 16S databases still label these as *Rhodococcus*.

## 3. Biosynthetic potential (antiSMASH 7.1.0, relaxed)
Files: [antiSMASH_BGC_summary.xlsx]({{artifact:874db76f-e23b-47bd-8a49-017504f23ea1}}), [antiSMASH_regions.csv]({{artifact:107585dd-e41f-4250-9007-7b68bb9b313f}}), and one full-result zip per genome.
- There are 57 regions in total: BF7 23, BF3 9, BF6 9, BF5 6, BF4 4, BF2 4, BF1 2.
- Only 3 regions reach ≥75% similarity to a known MIBiG cluster: ε-poly-L-lysine in BF7 and in BF4, and zeaxanthin in BF1. Most regions are below 50%, i.e. potentially novel.
- BF5 and BF4 (binned from the same culture) share NAPAA, T3PKS and betalactone region types.

![BGC types per genome]({{artifact:3216d9f7-d15c-401c-8371-7242e5a4bd11}})

## 4. Target genera in amplicon data
### 4a. Meta-analysis
2,411 samples with ≥1,000 reads. Methods: BLAST (≥97%) and exact V4 match against NCBI 16S references for 11 genera. Files: [amphibian_genus_abundance.xlsx]({{artifact:487b5319-6b32-4e83-9134-dc48b16f5f10}}), per-sample [amphibian_genus_relabund_by_sample.csv]({{artifact:7ce1bc36-0aaa-40ec-b461-d22e24b98de1}}).
- *Sphingomonas*: present in 80.9% of samples (mean 0.78%). *Arthrobacter* 41.7%, *Rhodococcus* 32.4%, *Kineococcus* 14.4%, *Ancylobacter* 13.7%.
- *Rhodococcoides* and *Psychromicrobium*: no sOTU has a BLAST top hit to either genus. The exact V4 matches that exist are shared with *Rhodococcus* or *Arthrobacter*.
- *Cryobacterium*: exact-V4 prevalence is 23.9%, falling to 16.9% when at least 2 supporting reference records are required.

### 4b. Your bullfrog dataset
42 samples, 3,700 ASVs of 120 bp; P10F3 and P4F1 were dropped for having <1,000 reads. Files: [bullfrog_genus_abundance.xlsx]({{artifact:9836893a-0a4b-41bb-ba47-0c6ae73c7590}}), [ASV_target_genus_matches.xlsx]({{artifact:a2be2d04-8839-402f-997d-b9252c2764f3}}).
- *Sphingomonas*: 0.93–1.42%, in 100% of samples. *Ancylobacter* about 0.08%.
- *Rhodococcoides* candidate (ASV 5cb28bd0; tied with *Antrihabitans*/*Rhodococcus*): 0.04%, in 50% of samples.
- *Kineococcus* 0.01–0.03%; *Arthrobacter* 0.008%.
- Not detected: *Psychromicrobium*. The *Cryobacterium* signal rests only on single-record matches and is likely spurious.

### 4c. Red-backed salamander (RBS1–15)
File: [salamander_genus_abundance.xlsx]({{artifact:dd16e100-c8c7-4cb2-b88e-cf5ea8204584}}).
- *Sphingomonas* 0.80% (part of this signal has BLAST top hits to *Parasphingomonas*). *Ancylobacter* 0.11%, in 53% of samples.
- Not detected: *Rhodococcoides*, *Psychromicrobium* and *Kineococcus*.

### 4d. Your isolates' own V4 sequences (from the genomes) in the meta-analysis
- BF5: 0.22% (rank 55). BF1: 0.052%. BF3: 0.012%. BF6: 0.0007%.
- Absent: BF7 and BF2.
- AmphiBac isolates that share BF7's V4 are mostly non-inhibitory (19 of 28). All 3 that share BF3's V4 are non-inhibitory.

## 5. AmphiBac 2025.1 (tested) vs the meta-analysis
Isolates were matched to sOTUs by exact containment of the 90-bp sOTU in the isolate 16S, otherwise by ≤3 mismatches over the V4 window (allowing for Sanger error). Files: [AmphiBac_meta_abundance.xlsx]({{artifact:33c8c7ee-0678-4428-92b2-fc3490b532e1}}), [AmphiBac_near_matches.xlsx]({{artifact:02bd023d-f792-4d52-84d4-34c94e8f3920}}).
- **Matching:** 6,610 of the 6,729 V4-spanning isolates match an sOTU with reads. 1,144 isolates lack V4 coverage (median length 488 bp).
- **Cultured share by abundance class:** 15/15 dominant, 55/82 intermediate, 494/48.3k rare sOTUs have a cultured representative. Per sample, cultured sOTUs are a median 20% of sOTUs but 68% of reads.
- **Detected vs undetected isolates:** inhibition rates are similar, 49% vs 43% (p = 0.32).

![Rank abundance]({{artifact:9b827439-f9ce-496c-a860-bf3339be6a37}})
![Cultured vs abundance]({{artifact:220050de-8253-43e7-8e00-2c1e2a620e6f}})

**Source host** ([AmphiBac_source_host_abundance.xlsx]({{artifact:1ce89640-3567-47af-87cb-80a6f015cea3}})):
- Isolates are usually less abundant on the host species they were cultured from than elsewhere: the source host is higher in only 26% of the 1,626 sOTU × host units (pooled).
- Host by host (1,476 units), the source host ranks above a median 54% of other host species.
- The taxon is absent on its source host in 35% of units, and the source host is the top host only 40 times.

![Source host]({{artifact:5e4f0031-2b62-4bca-861f-e61d7cbb5c07}})

## 6. Abundance class vs Bd function
Classes: dominant ≥1%, intermediate 0.1–1%, rare <0.1%. Files: [AmphiBac_abundance_class_function.xlsx]({{artifact:9fe895e2-f783-47fb-a6fe-01cb90cbfb07}}), [AmphiBac_abundance_class_definitions.xlsx]({{artifact:c6b6f3c7-566c-449f-bcb0-c9d33157accb}}).

| Definition | n (D / I / R) | % inhibitory (D / I / R) | p (χ²) |
|---|---|---|---|
| Global mean, taxa | 15 / 55 / 472 | 67 / 49 / 39 | 0.048 |
| Global mean, isolates | 1,017 / 1,797 / 3,158 | 68 / 52 / 40 | (inflated, non-independent) |
| Mean when present, taxa | 65 / 278 / 199 | 57 / 44 / 33 | 0.0013 |
| Host-specific (59 hosts; + 523 not detected) | 88 / 230 / 638 | 62 / 55 / 44 (not detected 42) | 0.0004 |

**Conclusion:** under every definition, rare taxa are *less* often Bd-inhibitory, not more. Framing for the paper: your isolates are an exception to the database-wide trend. They are rare or undetected yet strongly inhibitory, which shows that inhibition is a strain-level trait not predictable from 16S.

![Abundance class vs function]({{artifact:ba58bb37-446c-4c19-a534-3510f43ed1c0}})
![Alternative definitions]({{artifact:0f3418be-68ff-42f3-b136-11d8128ed69e}})

## 7. Variation in Bd inhibition within each sOTU
Score = (1 − proportional Bd growth) × 100, from 96-well cell-free supernatant (CFS) assays only. Current values, after restricting to the three Bd labels:
- **Spread:** 174 sOTUs have ≥5 scored isolates (4,551 isolates, 16 studies). The median within-sOTU SD is 39.6 points.
- **Opposite effects:** 74.7% of sOTUs contain both an isolate scoring >50 and one scoring <0 (enhancing Bd growth).
- **No class effect:** median SD by class is 37.0 / 38.8 / 40.6 (Kruskal–Wallis p = 0.92).
- **Label agreement** (sOTUs with ≥3 isolates): a median of 67% of isolates share the majority label, and 13% of sOTUs are unanimous.
- **Variance partitioning:** sOTU explains 0.15 of the variance (ICC), study 0.23, and sOTU × study 0.39. 30% of the 146 multi-study sOTUs show a significant study effect.

Current values from the rebuilt R project (script 05). Earlier saved versions of these three files predated the Bd-label filter and showed 175 sOTUs, 11% unanimous and 29% of 147; the current versions are linked here: [sOTU_Bd_function_variance.csv]({{artifact:96fbc9ef-037a-4586-b218-a6b87662d615}}), [Bd_variance_partitioning.csv]({{artifact:ea9d113e-c330-4e1f-b934-939d3c2c2e65}}), [sOTU_Bd_variance.png]({{artifact:107adaa7-705f-4d53-b289-4877bf7412f3}}).

**Panel A, packaged for the separate manuscript:** [sOTU_Bd_panelA.zip]({{artifact:c9243529-ebfa-4c5f-941f-c4ef3fb1cdb9}}) (data, R script, figure, methods; abundance colouring removed).
![Panel A]({{artifact:928aeecd-23f0-4a4b-a1d7-ee91c5fba0f9}})

## 8. Decisions recorded
- AmphiBac matching uses ≤3 mismatches over V4 (your choice, to allow for Sanger error).
- Taxon-level p-values are the conservative ones; isolate-level tests are inflated by non-independence.
- Genome and salamander analyses are kept out of the R package.

## 9. Multi-hit sensitivity of the abundance-class results
Files: [AmphiBac_multihit_sensitivity.xlsx]({{artifact:d8278fcb-ebcf-49e9-8385-ec27ffb7927f}}), figure [AmphiBac_multihit_sensitivity.png]({{artifact:32fcc817-b3d5-45af-ae89-10b90039c6fa}}); script code/multihit_sensitivity.py.
- **Isolate side:** 852 of the 5,972 matched Bd-tested isolates (14%) have no exact match and are tied between several nearest sOTUs (median 5). No isolate exactly contains more than one sOTU. The original "most abundant tied sOTU" rule placed 229 isolates in the dominant or intermediate class that a "rarest" rule would call rare.
- **sOTU side:** 75.5% of sOTUs with ≥2 isolates have conflicting Bd labels (dominant 13/13, intermediate 45/49, rare 189/265).
- **Taxon-level χ² (p = 0.048):** significant only under the most-abundant rule (0.048) or the all-tied rule (0.023). It is not significant with rarest-tied (0.26), drop-tied (0.29), exact-only (0.39) or random tie-breaking (median 0.16).
- **Direction:** dominant ≥ rare under every tie rule and under the majority, unanimous and random-isolate label rules. The exception is the "all isolates inhibitory" label rule, where rare > dominant (23.5% vs 6.7%, p = 0.008). That reversal reflects isolate count: an sOTU with 19 isolates rarely has all of them inhibitory.
- **Label rules:** "any isolate inhibitory" and "all isolates inhibitory" are confounded by isolate count (median 19 per dominant sOTU vs 2 per rare sOTU).
- **Robust to every rule:** isolate-level inhibitory share (68–69 / 51–52 / 40–41%), cultured share (dominant always 15/15), and a weak positive correlation between abundance and median score (ρ 0.12–0.14, p ≤ 0.009).
- **Manuscript wording:** "rare taxa are not more often inhibitory" stands; describe the dominant > rare gradient as a trend that is sensitive to how multi-hits are resolved.

## 10. Open items
1. **R project (`amphib_R`)**: rebuilt on 2026-11-02 by replaying the recorded session code, and now included in the download bundle. All outputs match the earlier results. The one exception is the *Arthrobacter* strictest-method count (unique to one genus, ≥2 references) in the bullfrog data: 0 ASVs vs 1 in the saved workbook, because the exact match is shared with *Paeniglutamicibacter*. The 16S BLAST was re-run against the 29 Oct 2026 NCBI release (27,646 sequences) with identical BLAST-based results.
2. The optional NCBI re-download script timed out against NCBI E-utilities and is untested. The reference sequences used are saved in [ref16S_11genera_ncbi.tar.gz]({{artifact:05d759da-dbb2-4d8a-84b1-fb5fc6866353}}).
3. Not re-checked for multi-hit sensitivity: the mean-when-present and host-specific class definitions and the source-host analysis, which use the same most-abundant tie rule.
4. Offered but not run: the ≥2-record *Cryobacterium* filter for the per-sample meta-analysis files, and bullfrog-only abundance classes.

![Multi-hit sensitivity]({{artifact:32fcc817-b3d5-45af-ae89-10b90039c6fa}})
