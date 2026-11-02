# session cell 101 | 2026-10-17 18:01 | env: python | language: python
# files written: []


S3=Sx|set(nd.best_sOTU[nd.min_mismatch<=3])
print("cultured sOTUs exact:",len(Sx),"<=3mm:",len(S3), "isolates rescued <=3:",int((nd.min_mismatch<=3).sum()))
rich_3=np.zeros(ns); rd_3=np.zeros(ns)
with open(fn) as f:
    f.readline()
    for line in f:
        fl=line.rstrip("\r\n").split("\t")
        if len(fl)!=ns+1 or len(fl[0])!=90 or fl[0] not in S3: continue
        v=np.array(fl[1:],dtype=np.float64); rich_3+=v>0; rd_3+=v
tmp=pd.DataFrame({"r3":rich_3,"d3":rd_3},index=samples).loc[gsamples]
PSm["cultured_sOTUs_3mm"]=tmp.r3; PSm["reads_cultured_3mm"]=tmp.d3
PSm["pct_sOTUs_cultured_3mm"]=100*PSm.cultured_sOTUs_3mm/PSm.richness; PSm["pct_reads_cultured_3mm"]=100*PSm.reads_cultured_3mm/PSm.total_reads
ST["cultured_3mm"]=ST.index.isin(S3); pres_=ST[ST.n_samples>0].copy(); pres_["cls"]=pres_.mean_rel_pct.apply(acls)
cf=pres_.groupby("cls").agg(n_sOTUs=("reads","size"),n_cultured_exact=("cultured_exact","sum"),n_cultured_1mm=("cultured_99","sum"),n_cultured_3mm=("cultured_3mm","sum"))
for c in ["exact","1mm","3mm"]: cf[f"pct_cultured_{c}"]=(100*cf[f"n_cultured_{c}"]/cf.n_sOTUs).round(2)
print(cf.to_string())
print(PSm[["pct_sOTUs_cultured_exact","pct_reads_cultured_exact","pct_sOTUs_cultured_3mm","pct_reads_cultured_3mm"]].median().round(1).to_dict())
d3=kn1.FastaID.isin(set(bt.FastaID[bt.detection.str.startswith("detected")])|set(nd.FastaID[nd.min_mismatch<=3]))
print("inhib frac detected(<=3)",round((kn1[d3].Bd_function=="Inhibitory").mean(),3),"not",round((kn1[~d3].Bd_function=="Inhibitory").mean(),3),"n not",int((~d3).sum()), sst.chi2_contingency(pd.crosstab(d3,kn1.Bd_function=="Inhibitory"))[1])
