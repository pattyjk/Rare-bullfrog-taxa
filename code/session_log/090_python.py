# session cell 90 | 2026-10-17 17:52 | env: python | language: python
# files written: []


presS=list(ST.index[ST.n_samples>0])
idx2=collections.defaultdict(list)
for s in presS: idx2[s[:K]].append(s)
def find2(seq):
    hits=set()
    for qq in (seq, seq.translate(rc)[::-1]):
        for i in range(len(qq)-89):
            lst=idx2.get(qq[i:i+K])
            if lst:
                w=qq[i:i+90]
                hits.update(s for s in lst if s==w)
    return hits
match2={k:find2(v) for k,v in fa.items()}
m["sOTUs_present"]=m.FastaID.map(lambda k:";".join(sorted(match2[k])))
m["has_V4"]=m.FastaID.map(lambda k:win[k] is not None)
m["detection"]=np.select([m.sOTUs_present!="",m.has_V4],["detected","not detected (V4 not in amphibian data)"],"unknown (no V4 coverage)")
print(m.detection.value_counts().to_dict())
bt=m[m.Bd_function.notna()]
print(pd.crosstab(bt.detection,bt.Bd_function,margins=True).to_string())
kn=bt[bt.detection!="unknown (no V4 coverage)"]
print("chi2 detected vs not x inhibitory:",sst.chi2_contingency(pd.crosstab(kn.detection,kn.Bd_function=="Inhibitory"))[:2], (kn.groupby("detection").Bd_function.apply(lambda x:(x=="Inhibitory").mean())).round(3).to_dict())
# sOTU-level redo
rows=[]
for _,x in bt[bt.sOTUs_present!=""].iterrows():
    for s in x.sOTUs_present.split(";"): rows.append((s,x.Bd_function,x.FastaID,x["Host.Species"],x.Country))
E=pd.DataFrame(rows,columns=["sOTU","Bd_function","FastaID","host","country"])
agg=E.groupby("sOTU").agg(n_isolates=("FastaID","nunique"),n_inhib=("Bd_function",lambda x:(x=="Inhibitory").sum()),n_non=("Bd_function",lambda x:(x=="Non-inhibitory").sum()),n_fac=("Bd_function",lambda x:(x=="Facilitating").sum()),
    n_hosts=("host","nunique"),countries=("country",lambda x:";".join(sorted(set(map(str,x))))))
agg["frac_inhib"]=agg.n_inhib/agg.n_isolates
agg["function_class"]=np.select([agg.frac_inhib>=0.5,agg.n_fac/agg.n_isolates>=0.5],["Inhibitory","Facilitating"],"Non-inhibitory/mixed")
agg=agg.join(ST[["reads","n_samples","prevalence_pct","mean_rel_pct","max_rel_pct"]])
agg["genus_meta"]=[tax.loc[s].genus_tax if s in tax.index else "" for s in agg.index]
agg["top_hit_meta"]=[str(tax.loc[s].Hit_name).split(" 16S")[0] if s in tax.index else "" for s in agg.index]
agg["abundance_rank"]=agg.index.map(rk); agg["percentile_rank"]=(100*agg.abundance_rank/len(present)).round(2)
agg["abundance_class"]=agg.mean_rel_pct.apply(acls)
print("sOTUs:",len(agg)); print(pd.crosstab(agg.abundance_class,agg.function_class,margins=True).to_string())
inh=agg[agg.function_class=="Inhibitory"]; non=agg[agg.function_class=="Non-inhibitory/mixed"]
t1=sst.mannwhitneyu(inh.mean_rel_pct,non.mean_rel_pct); t2=sst.mannwhitneyu(inh.prevalence_pct,non.prevalence_pct); t3=sst.spearmanr(agg.frac_inhib,agg.mean_rel_pct)
print(t1,t2,t3)
print(agg.groupby("function_class")[["mean_rel_pct","prevalence_pct","percentile_rank"]].median().round(4).to_string())
# rare vs dominant: % inhibitory sOTUs per class
print(agg.groupby("abundance_class").apply(lambda d:round(100*(d.function_class=="Inhibitory").mean(),1)).to_dict())
