# session cell 82 | 2026-10-17 17:49 | env: python | language: python
# files written: []


from scipy import stats as sst
mp=pd.read_csv(host.artifact_path("e39afe94-4d38-4b9a-bb95-7731a3f8da7e"),sep="\t",dtype=str).set_index("#SampleID")
gs=pd.Index(RM.columns)
# isolate -> sOTUs present
m["sOTUs_all"]=m.FastaID.map(lambda k:";".join(sorted(match[k])) if match[k] else "")
iso=m[m.Bd_function.notna()&(m.sOTUs_all!="")].copy()
print("Bd-tested isolates with sOTU match:",len(iso), iso.Bd_function.value_counts().to_dict())
# sOTU-level function aggregation
rows=[]
for _,x in iso.iterrows():
    for s in x.sOTUs_all.split(";"): rows.append((s,x.Bd_function,x.FastaID,x["Host.Species"],x.Country,x["Proportional.Growth.Bd"]))
E=pd.DataFrame(rows,columns=["sOTU","Bd_function","FastaID","host","country","pgrowth"])
agg=E.groupby("sOTU").agg(n_isolates=("FastaID","nunique"),n_inhib=("Bd_function",lambda x:(x=="Inhibitory").sum()),n_non=("Bd_function",lambda x:(x=="Non-inhibitory").sum()),n_fac=("Bd_function",lambda x:(x=="Facilitating").sum()),
    n_hosts=("host","nunique"),countries=("country",lambda x:";".join(sorted(set(map(str,x))))))
agg["frac_inhib"]=agg.n_inhib/agg.n_isolates
agg["function_class"]=np.select([agg.frac_inhib>=0.5,agg.n_fac/agg.n_isolates>=0.5],["Inhibitory (>=50% isolates)","Facilitating (>=50%)"],"Non-inhibitory/mixed")
agg=agg.join(ST[["reads","n_samples","prevalence_pct","mean_rel_pct","max_rel_pct"]],how="left")
agg[["reads","n_samples","prevalence_pct","mean_rel_pct","max_rel_pct"]]=agg[["reads","n_samples","prevalence_pct","mean_rel_pct","max_rel_pct"]].fillna(0)
agg["genus_meta"]=[tax.loc[s].genus_tax if s in tax.index else "" for s in agg.index]
agg["top_hit_meta"]=[str(tax.loc[s].Hit_name).split(" 16S")[0] if s in tax.index else "" for s in agg.index]
agg["pct_ident_meta"]=[tax.loc[s].Per_ident if s in tax.index else np.nan for s in agg.index]
# rank among all sOTUs present
present=ST[ST.n_samples>0]
rk=present.mean_rel_pct.rank(ascending=False,method="min")
agg["abundance_rank"]=agg.index.map(rk)
agg["percentile_rank"]=(100*agg.abundance_rank/len(present)).round(2)
def acls(x): return "absent" if x==0 else ("dominant (>=1%)" if x>=1 else ("intermediate (0.1-1%)" if x>=0.1 else "rare (<0.1%)"))
agg["abundance_class"]=agg.mean_rel_pct.apply(acls)
print("present sOTUs in data:",len(present))
ct=pd.crosstab(agg.abundance_class,agg.function_class,margins=True); print(ct.to_string())
a=agg[agg.mean_rel_pct>0]
for fc in a.function_class.unique(): print(fc, a[a.function_class==fc].mean_rel_pct.median().round(4), a[a.function_class==fc].prevalence_pct.median().round(2), (a.function_class==fc).sum())
inh=a[a.function_class.str.startswith("Inhib")]; non=a[a.function_class=="Non-inhibitory/mixed"]
print("MWU mean_rel inhib vs non:",sst.mannwhitneyu(inh.mean_rel_pct,non.mean_rel_pct))
print("MWU prevalence:",sst.mannwhitneyu(inh.prevalence_pct,non.prevalence_pct))
print("Spearman frac_inhib vs mean_rel:",sst.spearmanr(a.frac_inhib,a.mean_rel_pct))
# overall share of reads from AmphiBac-matched sOTUs
print("background: overall rare fraction of present sOTUs:", round(100*(present.mean_rel_pct<0.1).mean(),2))
