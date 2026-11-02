# session cell 129 | 2026-10-18 11:33 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_abundance_class_function.xlsx","sha256":"06d995532646a073d5ae2b00287ef574dc747bcac4c2ba4990f03e1699913519"}]


summ_i=pd.DataFrame({"abundance_class":order,"Bd_tested_isolates":n_i.values,"n_inhibitory":cti.Inhibitory.values,"n_no_effect":cti["No effect"].values,"n_enhancing":cti.Enhancing.values,
  "pct_inhibitory":pcti.Inhibitory.round(1).values,"pct_no_effect":pcti["No effect"].round(1).values,"pct_enhancing":pcti.Enhancing.round(1).values,
  "isolates_with_score":[len(d) for d in di],"median_score":[round(np.median(d),1) for d in di],"IQR_low":[round(np.percentile(d,25),1) for d in di],"IQR_high":[round(np.percentile(d,75),1) for d in di]})
tests_i=pd.DataFrame([dict(test="Isolate-level chi-square: function (3 categories) x abundance class",p=p_chi_i),dict(test="Isolate-level chi-square: inhibitory vs not x abundance class",p=p_chi2_i),
  dict(test="Isolate-level Kruskal-Wallis: inhibition score across classes",p=p_kw_i),
  dict(test="NOTE",p=np.nan,comment="Isolates of the same taxon/study are not independent; isolate-level p-values are inflated relative to the taxon-level tests.")])
Iout=I.merge(iso2[["FastaID","Isolate.ID","Host.Species","Country","Reference","Assay.Method"]],on="FastaID",how="left")[["FastaID","Isolate.ID","Host.Species","Country","Reference","Assay.Method","Bd_function","inhib_score","sOTU","a","abundance_class"]].rename(columns={"a":"sOTU_mean_rel_pct"})
with pd.ExcelWriter("AmphiBac_abundance_class_function.xlsx") as w:
    summ.to_excel(w,sheet_name="taxa_summary_by_class",index=False); tests.to_excel(w,sheet_name="taxa_tests",index=False)
    summ_i.to_excel(w,sheet_name="isolate_summary_by_class",index=False); tests_i.to_excel(w,sheet_name="isolate_tests",index=False)
    T_out.to_excel(w,sheet_name="taxa",index=False); Iout.to_excel(w,sheet_name="isolates",index=False); notes.to_excel(w,sheet_name="methods_notes",index=False)
print(summ_i.to_string())
