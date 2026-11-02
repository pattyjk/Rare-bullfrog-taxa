# session cell 31 | 2026-10-15 20:00 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphibian_genus_summary.csv","sha256":"40da9846e276b00a63f40d3da4ba5ed15593ac86d61402284e6b21021e156efe"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphibian_genus_abundance.xlsx","sha256":"86d8fa463c8e3ee7f76092be80ba64567bb7083da1663252bff9bd995b1dd6b6"}]


with pd.ExcelWriter("amphibian_genus_abundance.xlsx") as w:
    R.to_excel(w,sheet_name="genus_summary",index=False); S.to_excel(w,sheet_name="by_host_metadata",index=False)
    D.to_excel(w,sheet_name="sOTU_detail",index=False); notes.to_excel(w,sheet_name="methods_notes",index=False)
R.to_csv("amphibian_genus_summary.csv",index=False)
print(D[D.exact_V4_match_genera.str.contains("Rhodococcoides|Psychromicrobium")][["blast_top_hit","blast_pct_identity","exact_V4_match_genera","reads","n_samples","prevalence_pct","max_rel_pct"]].to_string())
print(R[R.method=="BLAST top hit, >=97% id"][["genus","n_sOTUs","reads","prevalence_pct","mean_rel_pct","median_rel_when_present_pct","max_rel_pct"]].to_string())
print(R[R.method.str.startswith("Exact V4 match to")][["genus","n_sOTUs","reads","prevalence_pct","mean_rel_pct","max_rel_pct"]].to_string())
