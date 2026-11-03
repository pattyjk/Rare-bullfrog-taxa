# Comparative genomics: BF7 vs *Rhodococcoides fascians* and BF1 vs *Sphingomonas paucimobilis*

Run date 2026-11-02. All tools capped at 4 threads (Panaroo ran with `-t 1`; see Caveats).

## 1. Reference genomes (`scripts/00_fetch_qc_ani_derep.py`)
- **Source.** NCBI Datasets REST API v2 (`api.ncbi.nlm.nih.gov/datasets/v2`), accessed 2026-11-02. The `datasets` CLI 18.40.0 failed TLS verification in the sandbox, so the REST API was used.
  - Reports were current assembly versions for taxid 1828 (*Rhodococcoides fascians*; the synonym *Rhodococcus fascians* resolves to the same taxid) and taxid 13689 (*S. paucimobilis*).
  - GCA/GCF pairs were collapsed to one row, preferring GCF: 214 R. fascians and 50 S. paucimobilis assemblies.
- **QC exclusions.**
  - NCBI atypical flag "contaminated".
  - CheckM completeness <90% where reported.
  - Additionally, NCBI atypical flag "genome length too small", because incomplete genomes inflate "absent in references" calls.
  - Excluded: 8 R. fascians and 16 S. paucimobilis assemblies. This includes the S. paucimobilis type strain assembly NCTC 11030 (GCF_900457515.1, CheckM 81.1%). The type strain is still represented by NBRC 13935 (GCF_000739895.2) and FDAARGOS_908 (GCF_016027095.1).
- **Species membership.** skani 0.4.0 `dist`, isolate vs. every QC-passing assembly. Members have ANI ≥95%; pairs skani does not report (<~80% ANI) were treated as non-members.
  - R. fascians: 193 of 206 are members. 7 have ANI 85.7–86.7% (NCBI best match *R. navarretei*); 6 were unreported.
  - S. paucimobilis: 26 of 34 are members; 8 were unreported (NCBI ANI check "Inconclusive", best matches in other genera or species).
- **Dereplication.** skani `triangle` among members. Greedy centroid clustering at ANI ≥99.9% with both aligned fractions ≥50%.
  - Representative priority: mandatory strains, then assembly level, then contig N50.
  - Single-linkage was tried first but chained 25 of 26 S. paucimobilis into one cluster, so it was rejected.
  - Mandatory strains were always kept: R. fascians DSM 43673 GCF_030813815.1 (MiGA closest), NBRC 12155 GCF_001894785.1 and LMG 3623 GCF_000761055.1 (type); S. paucimobilis FDAARGOS_881 GCF_016028435.1 (MiGA closest), NBRC 13935 GCF_000739895.2 and FDAARGOS_908 GCF_016027095.1 (type).
  - Result: 49 R. fascians and 10 S. paucimobilis references. The cap of 60 was not reached.
  - Every assembly, with status, metadata and ANI, is in `reference_genomes.tsv`. The accession lists are `accessions_used_*.txt`. Raw genome sets are not shipped.

## 2. Gene calling (`01_call_genes.py`)
- pyrodigal 3.7.1, single-genome training, translation table 11.
- The same procedure was applied to every genome, including BF7 and BF1 (their own assemblies `fg_inputs/<ID>/<ID>.fna`).
- Output was written as Prokka-style GFF3 with embedded FASTA for Panaroo.

## 3. Pangenomes (`03_pangenome_analysis.py`)
- **Primary method.** Panaroo 1.8.0 (CD-HIT 4.8.1), `--clean-mode strict --remove-invalid-genes`, default identity thresholds: `-c 0.98`, `--len_dif_percent 0.98`, family threshold 0.7.
- **Sensitivity method.** MMseqs2 18.8cc5c `easy-cluster --min-seq-id 0.5 -c 0.8 --cov-mode 0` (≥50% identity, ≥80% bidirectional coverage) on the same proteins.
- **Category thresholds.** Core ≥99% of genomes; soft core 95–99%; shell 15–95%; cloud <15%.
- **Isolate-specific definition.** Present in the isolate and in ≤10% of references: ≤4 of 49 for BF7, ≤1 of 10 for BF1.
- **Union and tiers.** The union of Panaroo-specific and MMseqs2-specific genes is reported with a `tier` column:
  - `both`: high confidence.
  - `mmseqs_only`: usually genes that Panaroo strict mode removed as rare or contig-edge content.
  - `panaroo_only`: the gene has a ≥50%-identity homolog in >10% of references, i.e. a divergent allele rather than absent content.
- **Low-coverage exclusion.** Genes on contigs with SPAdes coverage <5× were excluded and listed separately (`isolate_specific_lowcov_excluded_*.tsv`). In BF1 this removes 649 pyrodigal genes on ~1× low-GC fragments, including the Bacillus contig NODE_103 (3.4×).
- **Prokka mapping.** DIAMOND 2.2.8 blastp of specific pyrodigal proteins against the Prokka proteins in `fg_inputs`. Thresholds: ≥95% identity, ≥80% query and subject coverage, best hit. Prokka product, gene and EC were taken from the Prokka GFF.

## 4. Islands and domains (`05_pfam_islands.py`, `06_island_sharing.py`)
- **Pfam.** Pfam-A 38.2 (EBI current_release, downloaded 2026-11-02); HMMER 3.4 `hmmsearch --cut_ga`.
- **Islands.** Specific genes are grouped by contiguity on isolate contigs, allowing ≤2 intervening non-specific genes. An island has ≥3 specific genes; smaller segments are labelled `_seg###_small`.
- **BGC flags.** These are domain-content proxies only; antiSMASH was not run.
  - `likely_NRPS_BGC`: Condensation + AMP-binding domains in one island.
  - `likely_PKS_BGC`: ketoacyl-synthase (KS).
  - `likely_PKS/NRPS_hybrid`: KS + AMP-binding + PP-binding.
  - `possible_NRPS-like`: AMP-binding + PP-binding without a C domain.
  - Other classes flag single indicative domains: LANC_like, IucA_IucC, AfsA, terpene synthases.
  - KR/DH/AT domains alone were not used, because they also hit primary-metabolism enzymes.
- **Plasmid flag.** Putative plasmid contig = coverage ≥1.5× the length-weighted median coverage of contigs ≥10 kb (19.5× BF7, 28.0× BF1).
- **Island sharing.** For each island, the fraction of its genes whose MMseqs2 cluster contains each reference (`island_gene_sharing_matrix_*.tsv`). An island is called `novel_content` if ≥50% of its genes lack homologs (≥50% identity) in >10% of references.

## 5. Core-genome trees (`04_core_tree.py`, `07_tree_placement.py`)
- Single-copy core = Panaroo clusters present exactly once in every genome, excluding refound, pseudo and length-flagged genes. That gives 734 genes for R. fascians and 1,780 for S. paucimobilis.
- Per-gene nucleotide MAFFT 7.526 `--auto`, concatenation, then FastTree 2.2.0 `-nt -gtr -gamma`.
- Pairwise core SNPs were counted on ACGT-only sites. Trees are unrooted; placement used midpoint rooting. Labelled copies carry accession|strain|fas status.

## 6. fas / att / hyp virulence loci (`02_virulence_screen.py`)
- **Reference proteins.** From pFiD188 (GenBank JN093097.1; Francis et al. 2012):
  - fasR pFi_074, fasA–F pFi_077–pFi_82.
  - attR pFi_064, attX pFi_065, attA–H pFi_066–073.
  - hyp (putative RNA helicase) pFi_052.
- **Verification.** Gene assignments were confirmed by 100% identity and coverage of the original locus proteins Z29635.1 (fas operon) and AJ311775.1 (att locus).
- **Search.** DIAMOND blastp (e ≤1e-10) against pyrodigal proteins of BF7 plus all 193 conspecific assemblies (not only the 49 dereplicated).
  - A gene is present at ≥80% identity with query coverage ≥70% summed across ORFs; hyp is split into two ORFs in pyrodigal calls.
  - Observed identities were bimodal: ≥84.1% (orthologs) or ≤56.6% (chromosomal paralogs), with none in between.
- **Plasmid proxy.** The fraction of all 184 pFiD188 CDS present at the same thresholds is a proxy for a pFiD188-like plasmid.

## Caveats
- Presence of a gene or cluster means genetic capacity, not demonstrated expression or activity. BGC flags are domain-content proxies, not antiSMASH calls.
- Draft assemblies: BF7 has 219 contigs and BF1 has 2,456 contigs (142 ≥1 kb), and many references are drafts. Absence calls in references can reflect assembly gaps, and the R. fascians core (1,114 clusters at ≥99%) is depressed by fragmented references.
- Panaroo ran single-threaded because its multiprocessing pool raised a sandbox PermissionError. Results are unaffected; only the runtime (4 h for R. fascians) changed.
- Panaroo strict mode removes rare genes at contig edges, which under-calls isolate-specific content in fragmented assemblies. Hence the MMseqs2 union and the `tier` column. Filter to `tier == both` for the most conservative set.
- The S. paucimobilis reference set is small (10 after dereplication, 26 conspecific in total), so "≤10% of references" means ≤1 genome. BF1 is near-clonal with most references (44–85 core SNPs over 1.72 Mb), and its tree topology among those genomes is poorly supported.
- BF4/BF5 (mixed-culture bins) are not part of this track.
