# session cell 29 | 2026-10-15 20:00 | env: python | language: python
# files written: []


meta=m.loc[good]
strata=["host_order","lifestage","biogeo_realm","country","host_family"]
genmat={}
for lab,M in [("BLAST>=97",M1),("exactV4",M2)]:
    for g,ids in M.items():
        genmat[(lab,g)]=rel.loc[ids].sum() if ids else pd.Series(0.0,index=good)
long=[]
for st in strata:
    for lvl,idx in meta.groupby(st).groups.items():
        if len(idx)<20: continue
        for (lab,g),r in genmat.items():
            rr=r[idx]
            long.append(dict(stratum=st,level=lvl,n_samples=len(idx),method=lab,genus=g,prevalence_pct=round(100*(rr>0).mean(),2),mean_rel_pct=round(100*rr.mean(),4)))
S=pd.DataFrame(long)
for st in ["host_order","lifestage","biogeo_realm"]:
    x=S[(S.stratum==st)&(S.method=="BLAST>=97")].pivot(index="genus",columns="level",values="prevalence_pct")
    print(st, S[S.stratum==st].groupby("level").n_samples.first().to_dict()); print(x.to_string())
x=S[(S.stratum=="host_order")&(S.method=="exactV4")].pivot(index="genus",columns="level",values="prevalence_pct"); print(x.to_string())
