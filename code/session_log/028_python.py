# session cell 28 | 2026-10-15 20:00 | env: python | language: python
# files written: []


m=pd.read_csv("meta/map.txt",sep="\t",dtype=str).set_index("#SampleID")
T=totals.seq_reads
good=T.index[T>=1000]
print("samples used",len(good),"of",len(T))
rel=cnt[good].div(T[good],axis=1)
ti=tax.drop_duplicates("sOTU").set_index("sOTU")
present=cnt.index
# method 1: BLAST top hit (>=97% identity); method 1b: any identity
def m1(min_id):
    out={}
    for g in tgt6+["Rhodococcus (other)"]:
        ids=ti.index[(ti.target==g)&(ti.Per_ident>=min_id)]
        out[g]=[s for s in ids if s in present]
    return out
M1=m1(97); M1all=m1(0)
M2={g:[s for s in exm[exm.ref_genera.str.contains(g)].sOTU if s in present] for g in tgt6}
M2u={g:[s for s in exm[exm.ref_genera==g].sOTU if s in present] for g in tgt6}
def summ(ids):
    if not ids: return dict(n_sOTUs=0,reads=0,pct_all_reads=0,prevalence_pct=0,n_samples=0,mean_rel_pct=0,median_rel_when_present_pct=np.nan,max_rel_pct=0)
    r=rel.loc[ids].sum(); c=cnt.loc[ids,good].sum()
    p=r[r>0]
    return dict(n_sOTUs=len(ids),reads=int(c.sum()),pct_all_reads=round(100*c.sum()/T[good].sum(),4),
                prevalence_pct=round(100*(r>0).mean(),2),n_samples=int((r>0).sum()),mean_rel_pct=round(100*r.mean(),4),
                median_rel_when_present_pct=round(100*p.median(),4) if len(p) else np.nan,max_rel_pct=round(100*r.max(),2))
res=[]
for lab,M in [("BLAST top hit, >=97% id",M1),("BLAST top hit, any id",M1all),("Exact V4 match to NCBI 16S (any)",M2),("Exact V4 match, unique to genus",M2u)]:
    for g,ids in M.items():
        res.append(dict(genus=g,method=lab,**summ(ids)))
R=pd.DataFrame(res)
print(R.pivot(index="genus",columns="method",values="prevalence_pct").to_string())
print(R.pivot(index="genus",columns="method",values="mean_rel_pct").to_string())
print(R.pivot(index="genus",columns="method",values="n_sOTUs").to_string())
