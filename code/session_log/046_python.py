# session cell 46 | 2026-10-15 21:03 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/ASV_target_genus_matches.csv","sha256":"73fa1844f36c8ffb17cb30efe4b21a611b9971c03e73e00a55f15a3697a601b2"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/ASV_target_genus_matches.xlsx","sha256":"024037991e9a1ad49c2b65344a52c00b08b530252a2c31e3a1500bfc40e9cc9d"}]


# robust Cryobacterium amphib upper bound (>=2 supporting refs)
rob=[sid for sid in M2["Cryobacterium"] if sum(1 for v in refs["Cryobacterium"].values() if v==sid)>=2]
r=rel.loc[rob].sum(); print("Cryobacterium exactV4 >=2 refs: sOTUs",len(rob),"prevalence",round(100*(r>0).mean(),2),"mean",round(100*r.mean(),4))
# ASV exact120 support counts
C["exact120_n_refs"]=[sum(sum(c.values()) for c in ref120[s].values()) if s in ref120 else 0 for s in C.sequence]
C["exact120_single_ref_only"]=(C.exact120_n_refs==1)
D2=C[anyf].copy(); D2.index.name="ASV_ID"
cols2=["target_genera_any_method","top_hit_target_97","tied_targets","exact120_targets","exact120_n_refs","exact120_single_ref_only","blast_top_hit","blast_pct_identity","aln_length","tied_genera","n_tied_hits","tied_species","exact120_ref_species","in_meta_taxonomy_file","amphib_sOTU_match90","amphib_prevalence_pct","sequence"]
D2=D2[cols2].sort_values(["target_genera_any_method","blast_pct_identity"],ascending=[True,False])
allB=C[["blast_top_hit","blast_top_genus","blast_pct_identity","aln_length","tied_genera","n_tied_hits","exact120_ref_genera","amphib_sOTU_match90","sequence"]].copy()
nohit=pd.DataFrame({"blast_top_hit":"no hit (e<1e-10)","sequence":pd.Series(q)[[k for k in q if k not in C.index]]})
allB=pd.concat([allB,nohit]); allB.index.name="ASV_ID"
notes2=pd.DataFrame({"note":[
"Input: dna-sequences.fasta, 3,700 ASVs of 120 bp starting just after the 515F primer.",
"BLAST: blastn 2.17.0 -task blastn against NCBI 16S_ribosomal_RNA (27,647 sequences, downloaded 2026-10-16), evalue 1e-10, max 50 targets. 11 ASVs had no hit.",
"top_hit_target_97: genus of the first hit, with identity >=97% and alignment length >=100. tied_targets: target genera among all hits tied at the best bitscore (identity >=97%).",
"Rhodococcoides: a Rhodococcus or Rhodococcoides hit whose epithet is fascians, corynebacterioides, kroppenstedtii, kyotonensis/-e, navarretei, trifolii, yunnanensis/-e or luteus.",
"exact120: the 120 bp after 515F is identical to an NCBI nuccore 16S sequence (1000-1700 bp) of the genus. Downloaded for 10 genera only; Arthrobacter, Sphingomonas, Rhodococcus and Pseudarthrobacter capped at 1500 sequences.",
"exact120_single_ref_only = TRUE means only one GenBank record supports the match. These are often misidentified isolates (e.g. 'Cryobacterium sp.' matches whose BLAST top hit is 100% to a different genus's type strain), so treat them as unreliable.",
"amphib_sOTU_match90: the first 90 bp are identical to an sOTU in the amphibian meta-analysis table; amphib_prevalence_pct is that sOTU's prevalence across 2,411 amphibian samples."]})
with pd.ExcelWriter("ASV_target_genus_matches.xlsx") as w:
    SM.to_excel(w,sheet_name="summary",index=False); D2.to_excel(w,sheet_name="target_genus_ASVs"); allB.to_excel(w,sheet_name="all_ASV_blast"); notes2.to_excel(w,sheet_name="methods_notes",index=False)
D2.to_csv("ASV_target_genus_matches.csv")
print(D2.shape, allB.shape)
print(D2[["target_genera_any_method","exact120_n_refs","amphib_prevalence_pct"]].groupby("target_genera_any_method").agg(n=("exact120_n_refs","size"),single=("exact120_n_refs",lambda x:(x==1).sum()),in_amph=("amphib_prevalence_pct",lambda x:x.notna().sum())).to_string())
