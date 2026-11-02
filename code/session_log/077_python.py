# session cell 77 | 2026-10-17 16:58 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/salamander_genus_abundance.xlsx","sha256":"55b671dc1eccdeb5d5641cc6862382f7a9e5e93a8a5c93f546a2f64beaceff14"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/salamander_genus_relabund_by_sample.csv","sha256":"676c0a33d42900c7fffa85ff175a3fed831ccb73b0b5203dd5798071c16db258"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/salamander_genus_relabund_long.csv","sha256":"37800ee2947770b6d5882d05746844ad4829464515895bb302d0f39571398d43"}]


outW=pd.concat([pd.DataFrame({"total_reads":Tt.astype(int),"included_in_summary":inc}),W2],axis=1); outW.index.name="SampleID"
outW.to_csv("salamander_genus_relabund_by_sample.csv")
LB=W2.reset_index().rename(columns={"index":"SampleID"}).melt(id_vars="SampleID",var_name="gm",value_name="rel_abund_pct")
LB[["genus","method"]]=LB.gm.str.split("__",expand=True); LB=LB.drop(columns="gm"); LB["total_reads"]=LB.SampleID.map(Tt).astype(int)
LB.to_csv("salamander_genus_relabund_long.csv",index=False)
D["exact120_single_ref_only"]=D.exact120_n_refs==1
Dcols=["target_genera_any_method","methods","mean_rel_pct","prevalence_pct","max_rel_pct","blast_top_hit","blast_pct_identity","aln_length","tied_genera","n_tied_hits","tied_species","exact120_ref_genera","exact120_n_refs","exact120_single_ref_only","exact120_ref_species","metaTaxonomy90_hit","amphib_prevalence_pct","sequence"]
D=D[Dcols].sort_values(["target_genera_any_method","mean_rel_pct"],ascending=[True,False]); D.index.name="ASV_ID"
mdef=pd.DataFrame({"method":["BLAST97_top","BLAST97_tied","exact120","exact120_min2refs","exact120_unique_min2refs","metaTaxonomy90_97","any_method"],
 "definition":["The top BLAST hit against NCBI 16S RefSeq (27,647 seqs) is this genus (>=97% identity, alignment >=100 bp)",
  "The genus is among the hits tied for best bitscore (>=97% identity)",
  "The 120 bp after 515F is identical to an NCBI nuccore 16S sequence of the genus (11-genus reference set)",
  "exact120, supported by >=2 GenBank records",
  "exact120_min2refs, where the variant matches no other genus in the 11-genus reference set",
  "The first 90 bp are identical to an amphibian meta-analysis sOTU whose taxonomy is this genus at >=97%",
  "Union of BLAST97_top, BLAST97_tied, exact120 and metaTaxonomy90_97"]})
nts=pd.DataFrame({"note":["15 red-backed salamander samples (RBS1-15), 1,935 ASVs of 120 bp; 3,342-12,501 reads per sample. All samples included.",
 "Relative abundance = ASV reads / sample total x 100, summed per genus and method. Genus columns are not mutually exclusive.",
 "Rhodococcoides is assigned from the hit species epithet (fascians, corynebacterioides, kroppenstedtii, kyotonense, navarretei, trifolii, yunnanense, luteus). The reference set includes sequences filed under both old (Rhodococcus) and new (Rhodococcoides) names.",
 "Exact matches backed by a single GenBank record (exact120_single_ref_only=TRUE) are probably misidentified records: each such ASV is 100% identical to another genus's type strain (Staphylococcus, Glaciihabitans, Amnibacterium).",
 "The 'Arthrobacter' exact matches are to Arthrobacter sp./A. globiformis records, but BLAST gives 100% identity to Paeniglutamicibacter type strains, so they are ambiguous between the two genera."]})
with pd.ExcelWriter("salamander_genus_abundance.xlsx") as w:
    SB.to_excel(w,sheet_name="summary_long",index=False)
    SB.pivot(index="genus",columns="method",values="mean_rel_pct").to_excel(w,sheet_name="mean_rel_pct")
    SB.pivot(index="genus",columns="method",values="prevalence_pct").to_excel(w,sheet_name="prevalence_pct")
    outW.to_excel(w,sheet_name="per_sample_wide"); D.to_excel(w,sheet_name="ASVs")
    mdef.to_excel(w,sheet_name="method_definitions",index=False); nts.to_excel(w,sheet_name="notes",index=False)
print(D.shape, outW.shape, LB.shape)
print(SB[SB.method=="any_method"][["genus","n_ASVs","n_samples","median_rel_pct","max_rel_pct"]].to_string())
print(D[D.target_genera_any_method=="Sphingomonas"][["mean_rel_pct","prevalence_pct","blast_top_hit","blast_pct_identity"]].head(4).to_string())
