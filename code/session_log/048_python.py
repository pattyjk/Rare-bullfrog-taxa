# session cell 48 | 2026-10-15 21:10 | env: python | language: python
# files written: []


bf=pd.read_csv(host.artifact_path("b3169217-410e-4dad-bff1-4c925a06efb7"),sep="\t",index_col=0)
print(bf.shape, "ids in fasta:",bf.index.isin(list(q)).sum(), "fasta ids in table:",pd.Index(list(q)).isin(bf.index).sum())
Tb=bf.sum(); print(Tb.describe().round(0).to_dict(), (Tb<1000).sum())
RB=bf.div(Tb,axis=1)*100
ids_all=bf.index
Cx=C.reindex(ids_all)
# meta-analysis taxonomy via 90bp prefix
p90=pd.Series({k:q[k][:90] for k in ids_all})
def meta_gen(s,min_id=97):
    if s not in ti.index: return ""
    t=ti.loc[s]; return t.target if (isinstance(t.target,str) and t.Per_ident>=min_id) else ""
meta90=p90.apply(meta_gen)
def has(col,g): return Cx[col].fillna("").str.split(";").apply(lambda x:g in x)
M={}
for g in tg:
    blast97=(Cx.top_hit_target_97==g)
    tied=has("tied_targets",g)&(Cx.blast_pct_identity>=97)
    ex=has("exact120_targets",g)
    exr=ex&(Cx.exact120_n_refs>=2)
    exu=(Cx.exact120_targets==g)&(Cx.exact120_n_refs>=2)
    m90=(meta90==g)
    M[g]={"BLAST97_top":blast97,"BLAST97_tied":tied,"exact120":ex,"exact120_min2refs":exr,"exact120_unique_min2refs":exu,"metaTaxonomy90_97":m90,"any_method":blast97|tied|ex|m90}
cols={}; summ=[]
for g,d in M.items():
    gk=g.replace(" (other)","_other")
    for mk,mask in d.items():
        mask=mask.fillna(False).astype(bool)
        r=RB.loc[mask[mask].index].sum() if mask.any() else pd.Series(0.0,index=RB.columns)
        cols[f"{gk}__{mk}"]=r
        summ.append(dict(genus=g,method=mk,n_ASVs=int(mask.sum()),prevalence_pct=round(100*(r>0).mean(),1),n_samples=int((r>0).sum()),
                         mean_rel_pct=round(r.mean(),4),median_rel_pct=round(r.median(),4),max_rel_pct=round(r.max(),3)))
SB=pd.DataFrame(summ)
print(SB.pivot(index="genus",columns="method",values="mean_rel_pct").to_string())
print(SB.pivot(index="genus",columns="method",values="prevalence_pct").to_string())
print(SB.pivot(index="genus",columns="method",values="n_ASVs").to_string())
