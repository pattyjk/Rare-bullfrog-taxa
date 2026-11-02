# session cell 97 | 2026-10-17 17:59 | env: python | language: python
# files written: []


Sx=set(";".join(m.sOTUs_present[m.sOTUs_present!=""]).split(";"))
S1=Sx|set(nd.best_sOTU[nd.rescued_1mm])
print("cultured sOTUs exact:",len(Sx),"incl <=1mm:",len(S1))
rich=np.zeros(ns); rich_x=np.zeros(ns); rich_1=np.zeros(ns); rd_x=np.zeros(ns); rd_1=np.zeros(ns)
with open(fn) as f:
    f.readline()
    for line in f:
        fl=line.rstrip("\r\n").split("\t")
        if len(fl)!=ns+1 or len(fl[0])!=90: continue
        try: v=np.array(fl[1:],dtype=np.float64)
        except ValueError: continue
        p_=v>0; rich+=p_
        if fl[0] in Sx: rich_x+=p_; rd_x+=v
        if fl[0] in S1: rich_1+=p_; rd_1+=v
PSm=pd.DataFrame({"richness":rich,"cultured_sOTUs_exact":rich_x,"cultured_sOTUs_99":rich_1,"reads_cultured_exact":rd_x,"reads_cultured_99":rd_1,"total_reads":tot},index=samples).loc[gsamples]
PSm["pct_sOTUs_cultured_exact"]=100*PSm.cultured_sOTUs_exact/PSm.richness; PSm["pct_reads_cultured_exact"]=100*PSm.reads_cultured_exact/PSm.total_reads
PSm["pct_sOTUs_cultured_99"]=100*PSm.cultured_sOTUs_99/PSm.richness; PSm["pct_reads_cultured_99"]=100*PSm.reads_cultured_99/PSm.total_reads
PSm=PSm.join(mp[["host_order","host_species","lifestage","country"]])
print(PSm[["pct_sOTUs_cultured_exact","pct_reads_cultured_exact","pct_sOTUs_cultured_99","pct_reads_cultured_99"]].median().round(2).to_dict())
# class-level fraction cultured
ST["cultured_exact"]=ST.index.isin(Sx); ST["cultured_99"]=ST.index.isin(S1)
pres_=ST[ST.n_samples>0].copy(); pres_["cls"]=pres_.mean_rel_pct.apply(acls)
cf=pres_.groupby("cls").agg(n_sOTUs=("reads","size"),n_cultured_exact=("cultured_exact","sum"),n_cultured_99=("cultured_99","sum"))
cf["pct_cultured_exact"]=(100*cf.n_cultured_exact/cf.n_sOTUs).round(2); cf["pct_cultured_99"]=(100*cf.n_cultured_99/cf.n_sOTUs).round(2)
print(cf.to_string())
# finer bins for figure: by log10 mean abundance
