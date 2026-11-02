# session cell 117 | 2026-10-17 21:35 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_source_host_abundance.png","sha256":"f26e5c7c3650650edff631f323bb25b7d143530e87c9bcbcbd96e5b50824fe3f","preview_sha256":"881e061a26183720ea177131b377fcfaa23f896c610c0c8dd69696b683b649d0","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_005_AmphiBac_source_host_abundance.png","sha256":"585e2aff9af87577a906f8d6e084e4b6180fbc6d72d187674a19efb72a711b47","preview_sha256":"585e2aff9af87577a906f8d6e084e4b6180fbc6d72d187674a19efb72a711b47","preview_content_type":"image/png"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_source_host_abundance.xlsx","sha256":"e35a6d4d572be77fe478746f62d69e07fe5d1c62b954b6373aec59503a354b49"}]


a2.set_yticks(range(0,351,50))
fig.savefig("AmphiBac_source_host_abundance.png",dpi=300,bbox_inches="tight")
host.view_image("AmphiBac_source_host_abundance.png",crop=(0.5,0,1,0.3))
out=R2.drop(columns=["studies"]).copy()
out["studies"]=R2.studies
cols_o=["sOTUs","host_label","hsp","home_level","n_isolates","n_bd","n_inhib","function","home_n","home_mean","home_prev","other_mean","other_prev","log2_home_vs_other",
        "home_group_mean","median_other_groups","pct_other_groups_below_home","home_group_rank","n_groups","log2_home_vs_median_other_group","home_rank_among_species","n_host_species","top_host_species","studies"]
out=out[cols_o].rename(columns={"hsp":"host_species_normalised","home_mean":"mean_rel_pct_source_samples","other_mean":"mean_rel_pct_all_other_samples","home_prev":"prevalence_pct_source","other_prev":"prevalence_pct_other",
   "home_group_mean":"mean_rel_pct_source_group","median_other_groups":"median_of_other_group_means","n_isolates":"n_AmphiBac_isolates"})
summ=[]
for lv in ["same host species","same host genus","same host family"]:
    s=R2[R2.home_level==lv]; s2=s[s.home_group_mean.notna()]
    summ.append(dict(home_level=lv,n_units=len(s),pooled_median_source=round(s.home_mean.median(),4),pooled_median_other=round(s.other_mean.median(),4),pooled_pct_source_higher=round(100*(s.home_mean>s.other_mean).mean(),1),
        pooled_wilcoxon_p=sst.wilcoxon(s.home_mean,s.other_mean).pvalue,n_units_hostlevel=len(s2),hostlevel_pct_source_above_median_other=round(100*(s2.home_group_mean>s2.median_other_groups).mean(),1),
        hostlevel_median_pct_groups_below=round(s2.pct_other_groups_below_home.median(),1),hostlevel_pct_absent_on_source=round(100*(s2.home_group_mean==0).mean(),1),source_is_top_group=int((s2.home_group_rank==1).sum())))
SUMM=pd.DataFrame(summ)
fn_s=sp.groupby("function").agg(n=("function","size"),median_pct_below=("pct_other_groups_below_home","median"),pct_absent=("home_group_mean",lambda x:round(100*(x==0).mean(),1))).reset_index()
nts=pd.DataFrame({"note":[
"Isolates: AmphiBac 2025.1-tested, matched to amphibian meta-analysis sOTUs (exact containment, or <=3 mismatches over the 90-bp V4 window for Sanger error). Analysis unit = unique (matched sOTU set x source host species); isolates in a unit pooled, Bd function = Inhibitory if >=50% of Bd-tested isolates inhibitory.",
"Source host = AmphiBac Host.Species. If that species is not sampled in the meta-analysis (or is 'sp./aff.'), the closest relative is used: same host genus, else same host family (families for genera absent from the meta-analysis assigned manually; Gephryomantis/Ischconemia/Hypsiboas treated as Gephyromantis/Ischnocnema/Boana). 286 isolates with no relative in the meta-analysis excluded.",
"Pooled comparison: mean relative abundance in source-host samples vs mean across all other samples.",
"Host-level comparison: mean in the source host group vs the mean in each other host species/genus/family with >=5 samples; pct_other_groups_below_home = % of other host groups with strictly lower mean abundance.",
"2,411 samples with >=1,000 reads. Isolation and amplicon samples are generally from different studies/populations, so 'source host' is a host-taxon match, not the same animals."]})
with pd.ExcelWriter("AmphiBac_source_host_abundance.xlsx") as w:
    SUMM.to_excel(w,sheet_name="summary_by_level",index=False); fn_s.to_excel(w,sheet_name="species_level_by_function",index=False)
    out.to_excel(w,sheet_name="isolate_host_units",index=False); nts.to_excel(w,sheet_name="methods_notes",index=False)
print(SUMM.T.to_string()); print(fn_s.to_string())
