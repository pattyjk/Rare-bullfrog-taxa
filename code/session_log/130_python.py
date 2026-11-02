# session cell 130 | 2026-10-18 11:36 | env: python | language: python
# files written: []


assert {'E','iso2','so','nd','T','order'} <= globals().keys()
prev=so.set_index("sOTU").prevalence_pct.to_dict()
for s,p in zip(nd.best_sOTU,nd.best_sOTU_prevalence_pct): prev.setdefault(s,p)
T["prevalence_pct"]=T.index.map(prev)
T["mean_when_present_pct"]=T.mean_rel_pct/(T.prevalence_pct/100)
print(T[["mean_rel_pct","prevalence_pct","mean_when_present_pct"]].describe().round(3).to_string())
T["class_present"]=T.mean_when_present_pct.apply(acls)
print(pd.crosstab(T.abundance_class,T.class_present).reindex(index=order,columns=order).to_string())
ctp=pd.crosstab(T.class_present,T.function_class).reindex(order)[["Inhibitory","No effect / mixed","Enhancing"]]
print(ctp.to_string()); print((ctp.div(ctp.sum(axis=1),axis=0)*100).round(1).to_string())
dp=[T[(T.class_present==c)&T.median_inhib_score.notna()].median_inhib_score.values for c in order]
print([len(d) for d in dp],[round(np.median(d),1) for d in dp])
print("chi2 inhib vs not:",sst.chi2_contingency(np.c_[ctp.Inhibitory,ctp.sum(axis=1)-ctp.Inhibitory])[1],"3x3:",sst.chi2_contingency(ctp.values)[1],"KW:",sst.kruskal(*dp).pvalue)
print("Spearman mean_when_present vs frac inhib:",sst.spearmanr(T.mean_when_present_pct,T.n_inhib/T.n_isolates))
