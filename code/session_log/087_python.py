# session cell 87 | 2026-10-17 17:51 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_meta_abundance.xlsx","sha256":"470d140452127d01e215920930c4a8c637ffec26c7ba50fbebe816c53b58d06e"}]


taxa2={"BF7":"Rhodococcoides fascians","BF3":"Kineococcus sp. nov.","BF6":"Psychromicrobium sp. nov.","BF5":"Arthrobacter sp. nov.","BF2":"Xanthobacteraceae sp. nov.","BF1":"Sphingomonas paucimobilis"}
I2=ISO.copy()
I2.insert(1,"MiGA_taxonomy",I2.isolate.map(taxa2))
I2["note"]=""
I2.loc[(I2.isolate=="BF1")&(I2.V4_copy==2),"note"]="CONTAMINANT: Bacillus 16S on NODE_103 (3.4x coverage vs ~31x genome); exclude"
for lab_ in ["bullfrog","salamander"]:
    ex=I2[f"{lab_}_match"].str.startswith("exact")
    for c in ["prevalence_pct","mean_rel_pct","max_rel_pct"]:
        I2[f"{lab_}_nearestASV_{c}"]=np.where(ex,np.nan,I2[f"{lab_}_{c}"]); I2[f"{lab_}_{c}"]=np.where(ex,I2[f"{lab_}_{c}"],0.0)
I2["meta_abundance_class"]=I2.meta_mean_rel_pct.apply(acls)
cols_i=["isolate","MiGA_taxonomy","V4_copy","note","meta_taxonomy","meta_prevalence_pct","meta_mean_rel_pct","meta_max_rel_pct","meta_rank","meta_percentile","meta_abundance_class",
 "AmphiBac_isolates_same_V4","AmphiBac_Bd_function","AmphiBac_hosts","bullfrog_match","bullfrog_prevalence_pct","bullfrog_mean_rel_pct","bullfrog_max_rel_pct",
 "salamander_match","salamander_prevalence_pct","salamander_mean_rel_pct","salamander_max_rel_pct","bullfrog_nearestASV_prevalence_pct","bullfrog_nearestASV_mean_rel_pct","salamander_nearestASV_prevalence_pct","salamander_nearestASV_mean_rel_pct","V4_120bp"]
I2=I2[cols_i]
# summary tables
cross=pd.crosstab(agg.abundance_class,agg.function_class,margins=True)
fsum=a.groupby("function_class").agg(n_sOTUs=("n_isolates","size"),n_isolates=("n_isolates","sum"),median_mean_rel_pct=("mean_rel_pct","median"),median_prevalence_pct=("prevalence_pct","median"),median_percentile_rank=("percentile_rank","median")).reset_index()
tests=pd.DataFrame([dict(test="Mann-Whitney: mean rel. abundance, inhibitory vs non-inhibitory/mixed sOTUs",n=f"{len(inh)} vs {len(non)}",stat=float(sst.mannwhitneyu(inh.mean_rel_pct,non.mean_rel_pct).statistic),p=float(sst.mannwhitneyu(inh.mean_rel_pct,non.mean_rel_pct).pvalue)),
 dict(test="Mann-Whitney: prevalence, inhibitory vs non-inhibitory/mixed",n=f"{len(inh)} vs {len(non)}",stat=float(sst.mannwhitneyu(inh.prevalence_pct,non.prevalence_pct).statistic),p=float(sst.mannwhitneyu(inh.prevalence_pct,non.prevalence_pct).pvalue)),
 dict(test="Spearman: fraction of isolates inhibitory vs mean rel. abundance",n=str(len(a)),stat=float(sst.spearmanr(a.frac_inhib,a.mean_rel_pct).statistic),p=float(sst.spearmanr(a.frac_inhib,a.mean_rel_pct).pvalue))])
gsum=G.reset_index().rename(columns={"index":"sOTU"})[["sOTU","genus_meta","top_hit_meta","n_isolates","n_inhib","n_non","n_fac","frac_inhib","prevalence_pct","mean_rel_pct","max_rel_pct","percentile_rank","abundance_class"]].sort_values(["genus_meta","mean_rel_pct"],ascending=[True,False])
iso_out=m[["FastaID","Isolate.ID","Host.Species","Country","Life.Stage","Bd.Inhibition","Bd_function","Proportional.Growth.Bd","Reference","sOTUs_all"]].copy()
iso_out["n_sOTU_matches"]=iso_out.sOTUs_all.apply(lambda s:len(s.split(";")) if s else 0)
iso_out["sum_mean_rel_pct"]=iso_out.sOTUs_all.apply(lambda s:round(sum(ST.mean_rel_pct.get(x,0) for x in s.split(";")),5) if s else 0)
iso_out["max_prevalence_pct"]=iso_out.sOTUs_all.apply(lambda s:round(max(ST.prevalence_pct.get(x,0) for x in s.split(";")),2) if s else 0)
sotu_out=agg.reset_index().rename(columns={"index":"sOTU"}).sort_values("mean_rel_pct",ascending=False)
shr=share.reset_index().rename(columns={"index":"SampleID"})
notes4=pd.DataFrame({"note":[
"AmphiBac 2025.1 'tested' set: 7,873 isolates. Bd function normalised from Bd.Inhibition (Inhibitory/inhibitory/Inhibitor -> Inhibitory; NonInhibitory/NoEffect/No effect/Non-inhibitor/ns -> Non-inhibitory; Facilitating/enhancing -> Facilitating); 655 isolates with no Bd label excluded from function analyses.",
"Matching: an isolate matches an sOTU when the 90-bp sOTU occurs exactly within the isolate 16S sequence (either strand). 6,762 isolates matched 577 sOTUs. Some isolates match 2 sOTUs offset by ~19 bp (different study trimming); both are kept and summed.",
"Abundance: 2,411 amphibian meta-analysis samples with >=1,000 reads; relative abundance = sOTU reads / sample total. Classes by mean relative abundance across all samples: dominant >=1%, intermediate 0.1-1%, rare <0.1%, absent = 0 reads.",
"sOTU function_class: Inhibitory if >=50% of matched Bd-tested isolates inhibit Bd; Facilitating if >=50% facilitate; otherwise Non-inhibitory/mixed.",
"Your isolates: V4 (120 bp after 515F, followed by 806R site) extracted directly from the genome assemblies. BF4 contains no 16S gene. BF1 has a second, Bacillus-like 16S on a low-coverage contig (3.4x vs ~31x) that is a contaminant.",
"For bullfrog/salamander, abundance is reported only for an exact 120-bp ASV match; otherwise 0, with the nearest ASV (and its mismatches) given separately.",
"An exact V4 match means the same amplicon sequence, not necessarily the same strain."]})
with pd.ExcelWriter("AmphiBac_meta_abundance.xlsx") as w:
    I2.to_excel(w,sheet_name="your_isolates",index=False)
    fsum.to_excel(w,sheet_name="function_vs_abundance",index=False); cross.to_excel(w,sheet_name="class_crosstab"); tests.to_excel(w,sheet_name="tests",index=False)
    gsum.to_excel(w,sheet_name="target_genera_in_AmphiBac",index=False); sotu_out.to_excel(w,sheet_name="AmphiBac_sOTUs",index=False)
    iso_out.to_excel(w,sheet_name="AmphiBac_isolates",index=False); shr.to_excel(w,sheet_name="per_sample_share",index=False); notes4.to_excel(w,sheet_name="methods_notes",index=False)
print(fsum.to_string()); print(tests.to_string())
