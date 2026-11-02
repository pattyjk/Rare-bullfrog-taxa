# session cell 54 | 2026-10-15 21:23 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphibian_genus_relabund_by_sample.csv","sha256":"93b2fdc605f4b02f806fb41a43248f02d4be4bf8418802bedde71b082f2258c5"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphibian_genus_summary.csv","sha256":"79fdff1a44573f7cff22889e41cd4e6d52de60686f78add738e1cc5c50b4f85c"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphibian_genus_relabund_long.csv","sha256":"01ea43b29b71f27e899c6bfda332ef97fb3f4256010235185e1fdb11972c9961"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphibian_genus_abundance.xlsx","sha256":"ebba01ee50cdab9f756be03bfa21858acdbf1255c2c51d7cee757d25a6c06714"}]


def summ_amph(ids):
    ids=[s for s in ids if s in present]
    if not ids: return dict(n_sOTUs=0,reads=0,pct_all_reads=0,prevalence_pct=0,n_samples=0,mean_rel_pct=0,median_rel_when_present_pct=np.nan,max_rel_pct=0)
    r=rel.loc[ids].sum(); c=cnt.loc[ids,good].sum(); p=r[r>0]
    return dict(n_sOTUs=len(ids),reads=int(c.sum()),pct_all_reads=round(100*c.sum()/T[good].sum(),4),prevalence_pct=round(100*(r>0).mean(),2),n_samples=int((r>0).sum()),
                mean_rel_pct=round(100*r.mean(),4),median_rel_when_present_pct=round(100*p.median(),4) if len(p) else np.nan,max_rel_pct=round(100*r.max(),2))
R=R[R.genus!="Kineococcus"]
R=pd.concat([R,pd.DataFrame([dict(genus="Kineococcus",method=lab,**summ_amph(ids)) for lab,ids in Kmeth.items()])],ignore_index=True)
print(R[R.genus=="Kineococcus"].to_string())
# per-sample columns
relall=cnt.div(T.replace(0,np.nan),axis=1)*100
short={"BLAST top hit, >=97% id":"BLAST97","BLAST top hit, any id":"BLASTany","Exact V4 match to NCBI 16S (any)":"exactV4","Exact V4 match, unique to genus":"exactV4unique"}
for lab,ids in Kmeth.items():
    ids=[s for s in ids if s in present]
    out[f"Kineococcus__{short[lab]}"]=(relall.loc[ids].sum() if ids else 0.0).round(5)
out.to_csv("amphibian_genus_relabund_by_sample.csv")
Wk=out[[c for c in out.columns if "__" in c]]
L=Wk.reset_index().melt(id_vars="SampleID",var_name="gm",value_name="rel_abund_pct")
L[["genus","method"]]=L.gm.str.split("__",expand=True); L=L.drop(columns="gm"); L["total_reads"]=L.SampleID.map(T).astype(int)
L.to_csv("amphibian_genus_relabund_long.csv",index=False)
# strata
for lab,mk in [("BLAST>=97",k97),("exactV4",kex)]:
    ids=[s for s in mk if s in present]; r=rel.loc[ids].sum()
    for st in strata:
        for lvl,idx in meta.groupby(st).groups.items():
            if len(idx)<20: continue
            rr=r[idx]; S=pd.concat([S,pd.DataFrame([dict(stratum=st,level=lvl,n_samples=len(idx),method=lab,genus="Kineococcus",prevalence_pct=round(100*(rr>0).mean(),2),mean_rel_pct=round(100*rr.mean(),4))])],ignore_index=True)
# sOTU detail
kd=[]
for s in (kall|kex)&set(present):
    t=ti.loc[s]; r=rel.loc[s]; c=cnt.loc[s,good]
    kd.append(dict(sOTU=s,blast_top_hit=t.Hit_name,blast_pct_identity=t.Per_ident,blast_genus_assignment="Kineococcus" if s in kall else t.genus_tax,
        exact_V4_match_genera=";".join(sorted(v2g.get(s,[]))),reads=int(c.sum()),n_samples=int((r>0).sum()),prevalence_pct=round(100*(r>0).mean(),2),mean_rel_pct=round(100*r.mean(),5),max_rel_pct=round(100*r.max(),3)))
D=pd.concat([D[D.blast_genus_assignment!="Kineococcus"],pd.DataFrame(kd)],ignore_index=True)
with pd.ExcelWriter("amphibian_genus_abundance.xlsx") as w:
    R.to_excel(w,sheet_name="genus_summary",index=False); S.to_excel(w,sheet_name="by_host_metadata",index=False)
    D.to_excel(w,sheet_name="sOTU_detail",index=False); notes.to_excel(w,sheet_name="methods_notes",index=False)
R.to_csv("amphibian_genus_summary.csv",index=False)
print(S[(S.genus=="Kineococcus")&(S.method=="BLAST>=97")&S.stratum.isin(["host_order","lifestage","biogeo_realm"])][["stratum","level","prevalence_pct","mean_rel_pct"]].to_string())
# bullfrog Kineococcus ASVs in amphib
for a in kids: s=p90[a]; print(a[:8], s in present, round(100*(rel.loc[s]>0).mean(),2) if s in present else None)
print(out.shape)
