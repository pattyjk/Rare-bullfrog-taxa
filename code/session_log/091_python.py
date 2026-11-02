# session cell 91 | 2026-10-17 17:52 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_meta_abundance.xlsx","sha256":"4c033e91db0cf9e5c9ac384eca4a7155161ff820784cf349ef59c6fb61c692e5"}]


inh_s=agg.index[agg.function_class=="Inhibitory"]; non_s=agg.index[agg.function_class=="Non-inhibitory/mixed"]; fac_s=agg.index[agg.function_class=="Facilitating"]
need2=set(agg.index)-set(RM.index); print("sOTUs missing from RM:",len(need2))
RMx=RM.reindex(RM.index.union(list(need2)))
if need2:
    keep2={}
    with open(fn) as f:
        f.readline()
        for line in f:
            i_=line.find("\t")
            if line[:i_] in need2:
                fl=line.rstrip("\r\n").split("\t"); v=np.array(fl[1:],dtype=np.float64); keep2[fl[0]]=(v[good]/tot[good]*100).astype(np.float32)
    RMx=pd.concat([RM,pd.DataFrame.from_dict(keep2,orient="index",columns=gsamples)])
share=pd.DataFrame({"AmphiBac_matched_pct":RMx.loc[agg.index].sum(),"Bd_inhibitory_pct":RMx.loc[inh_s].sum(),"non_inhibitory_pct":RMx.loc[non_s].sum(),"facilitating_pct":RMx.loc[fac_s].sum()}).join(mp[["host_order","host_species","lifestage","country","biogeo_realm"]])
print(share.iloc[:,:4].median().round(2).to_dict(), share.groupby("host_order").Bd_inhibitory_pct.median().round(2).to_dict())
bfs=share[share.host_species.str.contains("catesbeian",na=False)]; print("bullfrogs:",len(bfs),bfs.iloc[:,:4].median().round(2).to_dict())
G=agg[agg.genus_meta.isin(tg7)]
gsum=G.reset_index().rename(columns={"index":"sOTU"})[["sOTU","genus_meta","top_hit_meta","n_isolates","n_inhib","n_non","n_fac","frac_inhib","prevalence_pct","mean_rel_pct","max_rel_pct","percentile_rank","abundance_class"]].sort_values(["genus_meta","mean_rel_pct"],ascending=[True,False])
print(G.groupby("genus_meta").agg(sOTUs=("n_isolates","size"),isolates=("n_isolates","sum"),inhib=("n_inhib","sum"),non=("n_non","sum"),fac=("n_fac","sum"),max_prev=("prevalence_pct","max")).to_string())
fsum=agg.groupby("function_class").agg(n_sOTUs=("n_isolates","size"),median_mean_rel_pct=("mean_rel_pct","median"),median_prevalence_pct=("prevalence_pct","median"),median_percentile_rank=("percentile_rank","median")).reset_index()
cls_pct=agg.groupby("abundance_class").apply(lambda d:pd.Series({"n_sOTUs":len(d),"n_inhibitory":int((d.function_class=="Inhibitory").sum()),"pct_inhibitory":round(100*(d.function_class=="Inhibitory").mean(),1)})).reset_index()
det=pd.crosstab(bt.detection,bt.Bd_function,margins=True)
c2=sst.chi2_contingency(pd.crosstab(kn.detection,kn.Bd_function=="Inhibitory"))
tests=pd.DataFrame([dict(test="Mann-Whitney, mean relative abundance: inhibitory vs non-inhibitory/mixed sOTUs",n=f"{len(inh)} vs {len(non)}",stat=float(t1.statistic),p=float(t1.pvalue)),
 dict(test="Mann-Whitney, prevalence: inhibitory vs non-inhibitory/mixed sOTUs",n=f"{len(inh)} vs {len(non)}",stat=float(t2.statistic),p=float(t2.pvalue)),
 dict(test="Spearman, fraction of isolates inhibitory vs sOTU mean relative abundance",n=str(len(agg)),stat=float(t3.statistic),p=float(t3.pvalue)),
 dict(test="Chi-square, isolate inhibitory (yes/no) vs detected / not detected in amphibian data",n=str(len(kn)),stat=float(c2[0]),p=float(c2[1]))])
iso_out=m[["FastaID","Isolate.ID","Host.Species","Country","Life.Stage","Bd.Inhibition","Bd_function","Proportional.Growth.Bd","Reference","detection","sOTUs_present"]].copy()
iso_out["sum_mean_rel_pct"]=iso_out.sOTUs_present.apply(lambda s:round(sum(ST.mean_rel_pct[x] for x in s.split(";")),5) if s else 0)
iso_out["max_prevalence_pct"]=iso_out.sOTUs_present.apply(lambda s:round(max(ST.prevalence_pct[x] for x in s.split(";")),2) if s else 0)
sotu_out=agg.reset_index().rename(columns={"index":"sOTU"}).sort_values("mean_rel_pct",ascending=False)
notes4.loc[1,"note"]="Matching: an isolate matches an sOTU when a 90-bp sOTU with reads in the amphibian table occurs exactly within the isolate 16S (either strand). Sequences in the taxonomy file with zero amphibian reads (incl. primer-containing/reverse-oriented artefacts) are excluded. detection: detected / not detected (isolate spans V4 but no amphibian sOTU matches) / unknown (isolate sequence does not cover V4)."
with pd.ExcelWriter("AmphiBac_meta_abundance.xlsx") as w:
    I2.to_excel(w,sheet_name="your_isolates",index=False)
    fsum.to_excel(w,sheet_name="function_vs_abundance",index=False); cls_pct.to_excel(w,sheet_name="pct_inhibitory_by_class",index=False)
    det.to_excel(w,sheet_name="detection_vs_function"); tests.to_excel(w,sheet_name="tests",index=False)
    gsum.to_excel(w,sheet_name="target_genera_in_AmphiBac",index=False); sotu_out.to_excel(w,sheet_name="AmphiBac_sOTUs",index=False)
    iso_out.to_excel(w,sheet_name="AmphiBac_isolates",index=False); share.reset_index().rename(columns={"index":"SampleID"}).to_excel(w,sheet_name="per_sample_share",index=False); notes4.to_excel(w,sheet_name="methods_notes",index=False)
print(cls_pct.to_string())
