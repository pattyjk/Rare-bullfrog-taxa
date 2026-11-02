# session cell 120 | 2026-10-18 11:29 | env: python | language: python
# files written: []


iso2=iso.merge(meta[["FastaID","Assay.Method"]],on="FastaID",how="left")
iso2["PG"]=pd.to_numeric(iso2["Proportional.Growth.Bd"],errors="coerce")
iso2["inhib_score"]=(1-iso2.PG)*100
ndm=nd[nd.min_mismatch<=3].set_index("FastaID").best_sOTU
iso2["sOTU_set"]=[ (s if isinstance(s,str) and s else (ndm.get(f,"") )) for f,s in zip(iso2.FastaID,iso2.sOTUs_present)]
abund=so.set_index("sOTU").mean_rel_pct.to_dict()
for s,a in zip(nd.best_sOTU,nd.best_sOTU_mean_rel_pct): abund.setdefault(s,a)
bt=iso2[iso2.Bd_function.notna()&(iso2.sOTU_set!="")]
rows=[(s,r.FastaID,r.Bd_function,r.inhib_score) for r in bt.itertuples() for s in r.sOTU_set.split(";")]
E=pd.DataFrame(rows,columns=["sOTU","FastaID","Bd_function","inhib_score"])
T=E.groupby("sOTU").agg(n_isolates=("FastaID","nunique"),n_inhib=("Bd_function",lambda x:(x=="Inhibitory").sum()),n_fac=("Bd_function",lambda x:(x=="Facilitating").sum()),
    n_non=("Bd_function",lambda x:(x=="Non-inhibitory").sum()),n_scored=("inhib_score",lambda x:x.notna().sum()),median_inhib_score=("inhib_score","median"),mean_inhib_score=("inhib_score","mean"))
T["mean_rel_pct"]=T.index.map(abund)
print("missing abundance:",T.mean_rel_pct.isna().sum())
def acls(x): return "Dominant (≥1%)" if x>=1 else ("Intermediate (0.1–1%)" if x>=0.1 else "Rare (<0.1%)")
T["abundance_class"]=T.mean_rel_pct.apply(acls)
T["function_class"]=np.select([T.n_inhib/T.n_isolates>=0.5,T.n_fac/T.n_isolates>=0.5],["Inhibitory","Enhancing"],"No effect / mixed")
order=["Dominant (≥1%)","Intermediate (0.1–1%)","Rare (<0.1%)"]
ct=pd.crosstab(T.abundance_class,T.function_class).reindex(order); print(ct.to_string())
print(T.groupby("abundance_class").median_inhib_score.describe().round(1).reindex(order).to_string())
sc=[T[(T.abundance_class==c)&T.median_inhib_score.notna()].median_inhib_score for c in order]
print("Kruskal:",sst.kruskal(*sc)); print("chi2 function x class:",sst.chi2_contingency(ct.fillna(0))[1], "fisher dominant vs rest inhibitory:")
# isolate-level scores by class (each isolate assigned to class of its most abundant sOTU)
iso_cls=E.assign(a=E.sOTU.map(abund)).sort_values("a",ascending=False).drop_duplicates("FastaID")
iso_cls["abundance_class"]=iso_cls.a.apply(acls)
print(iso_cls.groupby("abundance_class").inhib_score.agg(["count","median"]).reindex(order).round(1).to_string())
