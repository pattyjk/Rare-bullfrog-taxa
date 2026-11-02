# session cell 49 | 2026-10-15 21:11 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/bullfrog_genus_relabund_long.csv","sha256":"916508ac148f7f50ce32f118de7a95417dd4e8942cc80c41578900ffcbee104f"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/bullfrog_genus_relabund_by_sample.csv","sha256":"9c5f74b6467255ff7c35c9d57037dffaff3b66dc7545e73af3597ad6e52c3f08"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/bullfrog_genus_abundance.xlsx","sha256":"65ada8c9fc698127dedd6c135f1c768f20dab5e164c68fa9ca63ab039d1a50dd"}]


W2=pd.DataFrame(cols).round(5)
inc=Tb>=1000
summ=[]
for c in W2.columns:
    g,mk=c.split("__"); r=W2.loc[inc[inc].index,c]
    summ.append(dict(genus=g,method=mk,n_ASVs=int(M[g.replace("_other"," (other)")][mk].fillna(False).sum()),prevalence_pct=round(100*(r>0).mean(),1),n_samples=int((r>0).sum()),
                     mean_rel_pct=round(r.mean(),4),median_rel_pct=round(r.median(),4),max_rel_pct=round(r.max(),3)))
SB=pd.DataFrame(summ)
outW=pd.concat([pd.DataFrame({"total_reads":Tb.astype(int),"included_in_summary":inc}),W2],axis=1); outW.index.name="SampleID"
outW.to_csv("bullfrog_genus_relabund_by_sample.csv")
LB=W2.reset_index().rename(columns={"index":"SampleID"}).melt(id_vars="SampleID",var_name="gm",value_name="rel_abund_pct")
LB[["genus","method"]]=LB.gm.str.split("__",expand=True); LB=LB.drop(columns="gm"); LB["total_reads"]=LB.SampleID.map(Tb).astype(int)
# ASV-level table with per-method flags and abundance
fl=[]
for g,d in M.items():
    for mk,mask in d.items():
        for a in mask.fillna(False)[mask.fillna(False)].index: fl.append((a,g,mk))
F=pd.DataFrame(fl,columns=["ASV_ID","genus","method"])
asv=F.groupby(["ASV_ID","genus"]).method.apply(lambda x:";".join(x)).reset_index()
asv["mean_rel_pct"]=asv.ASV_ID.map(RB.loc[:,inc].mean(axis=1)).round(4); asv["prevalence_pct"]=asv.ASV_ID.map((RB.loc[:,inc]>0).mean(axis=1)*100).round(1)
asv["max_rel_pct"]=asv.ASV_ID.map(RB.loc[:,inc].max(axis=1)).round(3)
asv=asv.merge(Cx[["blast_top_hit","blast_pct_identity","tied_genera","exact120_ref_species","exact120_n_refs","amphib_prevalence_pct"]],left_on="ASV_ID",right_index=True,how="left")
asv["metaTaxonomy90_hit"]=asv.ASV_ID.map(lambda a:(lambda s:f"{ti.loc[s].Hit_name.split(' 16S')[0]} ({ti.loc[s].Per_ident}%)" if s in ti.index else "")(p90[a]))
asv["sequence"]=asv.ASV_ID.map(q)
asv=asv.sort_values(["genus","mean_rel_pct"],ascending=[True,False])
mnotes=pd.DataFrame({"method":["BLAST97_top","BLAST97_tied","exact120","exact120_min2refs","exact120_unique_min2refs","metaTaxonomy90_97","any_method"],
 "definition":["The top BLAST hit against NCBI 16S RefSeq is this genus (>=97% identity, alignment >=100 bp)",
  "The genus is among the hits tied for best bitscore (>=97% identity)",
  "The 120 bp after 515F is identical to any NCBI nuccore 16S sequence of the genus (10-genus reference set)",
  "exact120, supported by >=2 GenBank records (drops single, likely misidentified records)",
  "exact120_min2refs, where the variant matches no other genus in the 10-genus reference set",
  "The first 90 bp are identical to an amphibian meta-analysis sOTU whose meta-analysis taxonomy is this genus at >=97%",
  "Union of BLAST97_top, BLAST97_tied, exact120 and metaTaxonomy90_97 (most inclusive)"]})
gnotes=pd.DataFrame({"note":["42 samples; 2 samples with <1000 reads are flagged included_in_summary=FALSE and excluded from the summary sheet but kept in the per-sample file.",
 "Relative abundance = ASV reads / total sample reads x 100, summed over the ASVs assigned to the genus under each method.",
 "The one Rhodococcoides ASV (5cb28bd0...) is also identical to Rhodococcus sp. records, so its reads also appear under Rhodococcus_other for the exact120 methods. Genus columns are not mutually exclusive across methods.",
 "All Cryobacterium exact120 matches are supported by a single GenBank record each, and each of those ASVs is 100% identical to another genus's type strain, so they are probably spurious."]})
with pd.ExcelWriter("bullfrog_genus_abundance.xlsx") as w:
    SB.to_excel(w,sheet_name="summary_long",index=False)
    SB.pivot(index="genus",columns="method",values="mean_rel_pct").to_excel(w,sheet_name="mean_rel_pct")
    SB.pivot(index="genus",columns="method",values="prevalence_pct").to_excel(w,sheet_name="prevalence_pct")
    outW.to_excel(w,sheet_name="per_sample_wide"); asv.to_excel(w,sheet_name="ASVs",index=False)
    mnotes.to_excel(w,sheet_name="method_definitions",index=False); gnotes.to_excel(w,sheet_name="notes",index=False)
LB.to_csv("bullfrog_genus_relabund_long.csv",index=False)
print(Tb[~inc].to_dict())
print(SB[SB.method.isin(["BLAST97_top","BLAST97_tied","exact120_min2refs","any_method"])].pivot(index="genus",columns="method",values="mean_rel_pct").to_string())
print(SB[SB.method.isin(["BLAST97_top","BLAST97_tied","exact120_min2refs","any_method"])].pivot(index="genus",columns="method",values="prevalence_pct").to_string())
print(SB[SB.method=="any_method"][["genus","n_samples","median_rel_pct","max_rel_pct"]].to_string())
print(asv[asv.genus.isin(["Rhodococcoides","Arthrobacter"])][["genus","method","mean_rel_pct","prevalence_pct","max_rel_pct","blast_top_hit","metaTaxonomy90_hit"]].to_string())
