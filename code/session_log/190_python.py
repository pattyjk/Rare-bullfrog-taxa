# session cell 190 | 2026-11-02 11:16 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_multihit_sensitivity.xlsx","sha256":"f13454531f88b5431189c8881b71e1429ce1fb7e02df3f11941de4b5068c4cd5"}]

prev_tab.loc[prev_tab.quantity.str.contains('dominant / intermediate / rare') & prev_tab.quantity.str.contains('conflicting'),'value']='13/13 (100%) / 45/49 (92%) / 189/265 (71%)'
with pd.ExcelWriter('AmphiBac_multihit_sensitivity.xlsx') as w:
    for df,sh,ix in [(notes,'notes',False),(prev_tab,'multihit_prevalence',False),(R_,'taxon_pct_inhibitory',False),(rand_ci,'random_draw_ranges',True),
                     (SC,'taxon_score_by_tie_rule',False),(ISO,'isolate_level_by_tie_rule',False),(CUL,'cultured_share_by_tie_rule',False)]:
        df.to_excel(w,sheet_name=sh,index=ix)
print(prev_tab.to_string(index=False))