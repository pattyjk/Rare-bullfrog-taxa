# session cell 35 | 2026-10-15 20:55 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphibian_genus_relabund_by_sample.csv","sha256":"e3461839e881ded8c8b20ee9d4a77271aaf19cf7510122ccfc3a5bfa5125f003"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphibian_genus_relabund_long.csv","sha256":"b21bf411abef5b855e98c3ccfd64af0297ec22d0ee09e9fa9fe724497a517fb7"}]


T=totals.seq_reads
relall=cnt.div(T.replace(0,np.nan),axis=1)*100
methods={"BLAST97":M1,"BLASTany":M1all,"exactV4":M2,"exactV4unique":M2u}
cols={}
for mk,M in methods.items():
    for g,ids in M.items():
        gk=g.replace(" (other)","_other").replace(" ","_")
        cols[f"{gk}__{mk}"]=relall.loc[ids].sum() if ids else pd.Series(0.0,index=T.index)
W=pd.DataFrame(cols).round(5)
order=[f"{g}__{mk}" for g in ["Rhodococcoides","Rhodococcus_other","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter","Sphingomonas"] for mk in methods if f"{g}__{mk}" in W]
W=W[order]
mcols=["host_order","host_family","host_genus","host_species","lifestage","country","biogeo_realm","latitude_deg","longitude_deg","external_habitat_type","doi"]
out=pd.concat([pd.DataFrame({"total_reads":T.astype(int),"included_in_summary":T>=1000}),m.reindex(T.index)[mcols],W],axis=1)
out.index.name="SampleID"
out.to_csv("amphibian_genus_relabund_by_sample.csv")
# long format
L=W.reset_index().rename(columns={"index":"SampleID"}).melt(id_vars="SampleID",var_name="genus_method",value_name="rel_abund_pct")
L[["genus","method"]]=L.genus_method.str.split("__",expand=True); L=L.drop(columns="genus_method")
L["total_reads"]=L.SampleID.map(T).astype(int)
L[["SampleID","genus","method","rel_abund_pct","total_reads"]].to_csv("amphibian_genus_relabund_long.csv",index=False)
print(out.shape, L.shape, W.shape)
# consistency check vs earlier summary (mean over included samples)
chk=W.loc[T>=1000].mean().round(4)
print(chk[["Sphingomonas__BLAST97","Arthrobacter__BLAST97","Rhodococcoides__exactV4","Cryobacterium__exactV4"]].to_dict())
print("NaN samples (0 reads):",int(W.isna().any(axis=1).sum()))
