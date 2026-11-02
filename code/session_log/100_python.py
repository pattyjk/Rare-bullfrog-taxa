# session cell 100 | 2026-10-17 18:00 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_cultured_vs_abundance.png","sha256":"ef2335c51ebee749c45cb57b47228c00fc163d3d49d89666b34883d273279798","preview_sha256":"3f207a68abbbbc56928ba08daf1db7abdbe9bc58974a6778a793b7635750b03c","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_near_matches.xlsx","sha256":"97bb56884c1e8de16da6f002e9d862a1e7d447cdbd1f34c9956e1f94300076ab"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_005_AmphiBac_cultured_vs_abundance.png","sha256":"7b4d19b2c5fecdbd2fa633f1db9c2043969b6f4d0866b8fe747937593812f9a4","preview_sha256":"7b4d19b2c5fecdbd2fa633f1db9c2043969b6f4d0866b8fe747937593812f9a4","preview_content_type":"image/png"}]


for t in list(axs[0].texts):
    if t.get_text().startswith(("428","476")): t.remove()
axs[0].text(2,4,"428/48.3k (0.9%)\n476 (1.0%)",ha="center",va="bottom",fontsize=6)
fig.savefig("AmphiBac_cultured_vs_abundance.png",dpi=300,bbox_inches="tight")
host.view_image("AmphiBac_cultured_vs_abundance.png",crop=(0,0,0.5,1))
# save near-match table into workbook additions
nd_out=nd[["FastaID","Isolate.ID","Host.Species","Country","Bd_function","Reference","min_mismatch","pct_identity","n_tied_sOTUs","best_sOTU"]].copy()
nd_out["best_sOTU_taxonomy"]=nd_out.best_sOTU.map(lambda s:str(tax.loc[s].Hit_name).split(" 16S")[0] if s in tax.index else "")
nd_out["best_sOTU_mean_rel_pct"]=nd_out.best_sOTU.map(ST.mean_rel_pct).round(5); nd_out["best_sOTU_prevalence_pct"]=nd_out.best_sOTU.map(ST.prevalence_pct).round(2)
mmsum=pd.DataFrame({"max_mismatches":[0,1,2,3,5],"identity_threshold":["100%","≥98.9%","≥97.8%","≥96.7%","≥94.4%"],
   "isolates_detected":[5392]+[5392+cum[t] for t in [1,2,3,5]],"of_V4_spanning_isolates":6729})
mmsum["pct_detected"]=(100*mmsum.isolates_detected/mmsum.of_V4_spanning_isolates).round(1)
with pd.ExcelWriter("AmphiBac_near_matches.xlsx") as w:
    mmsum.to_excel(w,sheet_name="threshold_summary",index=False); nd_out.to_excel(w,sheet_name="unmatched_V4_isolates",index=False)
    cf.reset_index().rename(columns={"cls":"abundance_class"}).to_excel(w,sheet_name="cultured_by_class",index=False)
    PSm.reset_index().rename(columns={"index":"SampleID"}).to_excel(w,sheet_name="per_sample_cultured_share",index=False)
print(mmsum.to_string())
