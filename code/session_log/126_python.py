# session cell 126 | 2026-10-18 11:32 | env: python | language: python
# files written: []


assert {'E','abund','order','iso_cls'} <= globals().keys()
I=iso_cls.copy()
I["func"]=I.Bd_function.replace({"Non-inhibitory":"No effect","Facilitating":"Enhancing"})
n_i=I.abundance_class.value_counts().reindex(order)
cti=pd.crosstab(I.abundance_class,I.func).reindex(order)[["Inhibitory","No effect","Enhancing"]]
pcti=cti.div(cti.sum(axis=1),axis=0)*100
print(n_i.to_dict()); print(cti.to_string()); print(pcti.round(1).to_string())
di=[I[(I.abundance_class==c)&I.inhib_score.notna()].inhib_score.values for c in order]
print([len(d) for d in di],[round(np.median(d),1) for d in di])
p_chi_i=sst.chi2_contingency(cti.values)[1]; p_chi2_i=sst.chi2_contingency(np.c_[cti.Inhibitory,cti.sum(axis=1)-cti.Inhibitory])[1]; p_kw_i=sst.kruskal(*di).pvalue
print(p_chi_i,p_chi2_i,p_kw_i)
print("isolates mapping to >1 sOTU:", int(E.groupby("FastaID").sOTU.nunique().gt(1).sum()))
