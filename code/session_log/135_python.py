# session cell 135 | 2026-10-18 11:38 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_abundance_class_definitions.xlsx","sha256":"c22e794f2a4273f8545503219cedf2ce98d9e194b46a58228499e8847a690ffe"}]


SR=pd.DataFrame(stats_rows)
cross=pd.crosstab(T.abundance_class,T.class_present).reindex(index=order,columns=order)
cross2=pd.crosstab(U.global_class,U.host_class).reindex(index=order,columns=order4)
Tx=T.reset_index()[["sOTU","n_isolates","n_inhib","n_non","n_fac","function_class","median_inhib_score","mean_rel_pct","prevalence_pct","mean_when_present_pct","abundance_class","class_present"]].rename(columns={"abundance_class":"class_global","class_present":"class_mean_when_present"})
Ux=U.rename(columns={"hsp":"source_host_species","host_class":"class_on_source_host","global_class":"class_global"})
hsum=pd.DataFrame({"class":order4,"units":cth.sum(axis=1).values,"pct_inhibitory":(cth.Inhibitory/cth.sum(axis=1)*100).round(1).values,
   "isolates":cthi.sum(axis=1).values,"isolate_pct_inhibitory":(cthi.Inhibitory/cthi.sum(axis=1)*100).round(1).values,
   "isolate_median_score":Iu.groupby("host_class").inhib_score.median().reindex(order4).round(1).values})
nts=pd.DataFrame({"note":[
 "All three definitions use the same 0.1% / 1% thresholds; only the abundance statistic differs.",
 "Global mean: sOTU mean relative abundance across all 2,411 meta-analysis samples, zeros included (original definition).",
 "Mean when present: mean relative abundance across only the samples in which the sOTU is detected (= global mean / prevalence). Separates 'rare everywhere it occurs' from 'patchy but abundant where present'.",
 "Host-specific: for each isolate, the sOTU's mean relative abundance (zeros included) across meta-analysis samples of the isolate's source host species. Only source host species with >=5 meta-analysis samples (59 species; 3,253 isolates). Unit = sOTU x source host species. 'Not detected' = zero reads in every sample of that host species.",
 "Function per unit: Inhibitory if >=50% of Bd-tested isolates inhibitory, Enhancing if >=50% facilitating, else No effect/mixed. Inhibition score = (1 - proportional Bd growth) x 100 from 96-well CFS assays; unit value = median of its isolates.",
 "chi-square tests compare inhibitory vs not across the three detected classes (Not detected excluded); Kruskal-Wallis across the three detected classes."]})
with pd.ExcelWriter("AmphiBac_abundance_class_definitions.xlsx") as w:
    SR.to_excel(w,sheet_name="summary_all_definitions",index=False); hsum.to_excel(w,sheet_name="host_specific_summary",index=False)
    cross.to_excel(w,sheet_name="global_vs_when_present"); cross2.to_excel(w,sheet_name="global_vs_host_specific")
    Tx.to_excel(w,sheet_name="taxa_global_when_present",index=False); Ux.to_excel(w,sheet_name="units_host_specific",index=False); nts.to_excel(w,sheet_name="methods_notes",index=False)
print(SR.to_string())
