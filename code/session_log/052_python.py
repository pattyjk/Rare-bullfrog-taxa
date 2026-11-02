# session cell 52 | 2026-10-15 21:23 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/bullfrog_genus_relabund_long.csv","sha256":"dddc70f8d2c1964b5b7c4e80a5b391e50bd73233533d52779dbb672421bf692e"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/bullfrog_genus_relabund_by_sample.csv","sha256":"e65d66b6256d53957700da5c57538010ff82af6a7a3041411020ea69a965d3ce"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/bullfrog_genus_abundance.xlsx","sha256":"43e9571709b16db79a29d8cd5199ad997aaee2ff3a6eadf9d8ff13b5990ac982"}]


for mk,mask in M["Kineococcus"].items():
    mask=mask.fillna(False).astype(bool)
    W2[f"Kineococcus__{mk}"]=(RB.loc[mask[mask].index].sum() if mask.any() else pd.Series(0.0,index=RB.columns)).round(5)
summ=[]
for c in W2.columns:
    gname,mk=c.split("__"); r=W2.loc[inc[inc].index,c]
    summ.append(dict(genus=gname,method=mk,n_ASVs=int(M[gname.replace("_other"," (other)")][mk].fillna(False).sum()),prevalence_pct=round(100*(r>0).mean(),1),n_samples=int((r>0).sum()),
                     mean_rel_pct=round(r.mean(),4),median_rel_pct=round(r.median(),4),max_rel_pct=round(r.max(),3)))
SB=pd.DataFrame(summ)
outW=pd.concat([pd.DataFrame({"total_reads":Tb.astype(int),"included_in_summary":inc}),W2],axis=1); outW.index.name="SampleID"
outW.to_csv("bullfrog_genus_relabund_by_sample.csv")
LB=W2.reset_index().rename(columns={"index":"SampleID"}).melt(id_vars="SampleID",var_name="gm",value_name="rel_abund_pct")
LB[["genus","method"]]=LB.gm.str.split("__",expand=True); LB=LB.drop(columns="gm"); LB["total_reads"]=LB.SampleID.map(Tb).astype(int)
LB.to_csv("bullfrog_genus_relabund_long.csv",index=False)
# ASV rows for Kineococcus
kr=[]
for a in kids:
    meths=";".join(mk for mk,mask in M["Kineococcus"].items() if bool(mask.get(a,False)))
    kr.append(dict(ASV_ID=a,genus="Kineococcus",method=meths))
kr=pd.DataFrame(kr)
kr["mean_rel_pct"]=kr.ASV_ID.map(RB.loc[:,inc].mean(axis=1)).round(4); kr["prevalence_pct"]=kr.ASV_ID.map((RB.loc[:,inc]>0).mean(axis=1)*100).round(1)
kr["max_rel_pct"]=kr.ASV_ID.map(RB.loc[:,inc].max(axis=1)).round(3)
kr=kr.merge(B[["blast_top_hit","blast_pct_identity","tied_genera"]],left_on="ASV_ID",right_index=True,how="left")
kr["exact120_n_refs"]=0; kr["exact120_ref_species"]=""
kr["amphib_prevalence_pct"]=np.nan
kr["metaTaxonomy90_hit"]=kr.ASV_ID.map(lambda a:(lambda s:f"{ti.loc[s].Hit_name.split(' 16S')[0]} ({ti.loc[s].Per_ident}%)" if s in ti.index else "")(p90[a]))
kr["sequence"]=kr.ASV_ID.map(q)
asv=pd.concat([asv[asv.genus!="Kineococcus"],kr[asv.columns]]).sort_values(["genus","mean_rel_pct"],ascending=[True,False])
with pd.ExcelWriter("bullfrog_genus_abundance.xlsx") as w:
    SB.to_excel(w,sheet_name="summary_long",index=False)
    SB.pivot(index="genus",columns="method",values="mean_rel_pct").to_excel(w,sheet_name="mean_rel_pct")
    SB.pivot(index="genus",columns="method",values="prevalence_pct").to_excel(w,sheet_name="prevalence_pct")
    outW.to_excel(w,sheet_name="per_sample_wide"); asv.to_excel(w,sheet_name="ASVs",index=False)
    mnotes.to_excel(w,sheet_name="method_definitions",index=False); gnotes.to_excel(w,sheet_name="notes",index=False)
print(SB[SB.genus=="Kineococcus"].to_string())
print(kr[["ASV_ID","method","mean_rel_pct","prevalence_pct","max_rel_pct","blast_top_hit","blast_pct_identity","metaTaxonomy90_hit"]].to_string())
