# session cell 132 | 2026-10-18 11:37 | env: python | language: python
# files written: []


import re
mp=pd.read_csv(host.artifact_path("e39afe94-4d38-4b9a-bb95-7731a3f8da7e"),sep="\t",dtype=str).set_index("#SampleID").loc[RMh.columns]
mp["hsp"]=mp.host_genus.str.lower().str.strip()+" "+mp.host_species.str.lower().str.strip()
spN=mp.hsp.value_counts(); ok=set(spN[spN>=5].index)-{h for h in spN.index if h.endswith((" sp"," aff"," cf"))}
fix={"gephryomantis":"gephyromantis","ischconemia":"ischnocnema","hypsiboas":"boana"}
def nh(s):
    w=re.sub(r"\s+"," ",str(s).replace("."," ").replace("_"," ").strip().lower()).split()[:2]
    if w: w[0]=fix.get(w[0],w[0])
    return " ".join(w)
hostmap=iso2.set_index("FastaID")["Host.Species"].map(nh)
E2=E.copy(); E2["hsp"]=E2.FastaID.map(hostmap)
E2=E2[E2.hsp.isin(ok)]
print("isolates with source host sampled (>=5 samples):",E2.FastaID.nunique(),"host species:",E2.hsp.nunique())
# host-specific mean per (sOTU, host)
hm=RMh.T.groupby(mp.hsp).mean().T   # sOTU x host species
hp=(RMh.T>0).groupby(mp.hsp).mean().T*100
U=E2.groupby(["sOTU","hsp"]).agg(n_isolates=("FastaID","nunique"),n_inhib=("Bd_function",lambda x:(x=="Inhibitory").sum()),n_fac=("Bd_function",lambda x:(x=="Facilitating").sum()),
    median_inhib_score=("inhib_score","median")).reset_index()
U["host_mean_pct"]=[hm.at[s,h] for s,h in zip(U.sOTU,U.hsp)]; U["host_prev_pct"]=[hp.at[s,h] for s,h in zip(U.sOTU,U.hsp)]
U["host_n_samples"]=U.hsp.map(spN)
order4=["Dominant (≥1%)","Intermediate (0.1–1%)","Rare (<0.1%)","Not detected on source host"]
U["host_class"]=np.where(U.host_mean_pct==0,"Not detected on source host",U.host_mean_pct.apply(acls))
U["function_class"]=np.select([U.n_inhib/U.n_isolates>=0.5,U.n_fac/U.n_isolates>=0.5],["Inhibitory","Enhancing"],"No effect / mixed")
U["global_class"]=U.sOTU.map(T.abundance_class)
print(pd.crosstab(U.global_class,U.host_class).reindex(index=order,columns=order4).to_string())
cth=pd.crosstab(U.host_class,U.function_class).reindex(order4)[["Inhibitory","No effect / mixed","Enhancing"]]
print(cth.to_string()); print((cth.div(cth.sum(axis=1),axis=0)*100).round(1).to_string())
dh=[U[(U.host_class==c)&U.median_inhib_score.notna()].median_inhib_score.values for c in order4]
print([len(d) for d in dh],[round(np.median(d),1) if len(d) else None for d in dh])
c3=cth.loc[order]; print("chi2 inhib vs not (detected classes):",sst.chi2_contingency(np.c_[c3.Inhibitory,c3.sum(axis=1)-c3.Inhibitory])[1],"incl not-detected:",sst.chi2_contingency(np.c_[cth.Inhibitory,cth.sum(axis=1)-cth.Inhibitory])[1])
print("KW (4 groups):",sst.kruskal(*[d for d in dh if len(d)]).pvalue, "KW detected only:",sst.kruskal(*dh[:3]).pvalue)
# isolate-level
Iu=E2.merge(U[["sOTU","hsp","host_class"]],on=["sOTU","hsp"])
cthi=pd.crosstab(Iu.host_class,Iu.Bd_function).reindex(order4); print(cthi.to_string()); print((cthi.div(cthi.sum(axis=1),axis=0)*100).round(1).to_string())
print(Iu.groupby("host_class").inhib_score.agg(["count","median"]).reindex(order4).round(1).to_string())
