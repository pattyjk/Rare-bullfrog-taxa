# session cell 122 | 2026-10-18 11:30 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_002_AmphiBac_abundance_class_function.png","sha256":"6c426fc978242ea3b6ca2b1642eee4deb5e2c799fa8bd2521c0fd7b60ce0b7f0","preview_sha256":"6c426fc978242ea3b6ca2b1642eee4deb5e2c799fa8bd2521c0fd7b60ce0b7f0","preview_content_type":"image/png"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_abundance_class_function.xlsx","sha256":"2a7083af98c35a2ea33cff7cc4f36344ce3838fab2d15a7183ea8e0d940d132a"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_abundance_class_function.png","sha256":"235e37a3cc759e0fe33183c1ba88e97a4ce791a6e6e1d93335ffeb5a317b850e","preview_sha256":"9be113b0d67722f98a45848e026c08b944ae175a436b44fd53557ec46e04d546","preview_content_type":"image/jpeg"}]


inh_tab=pd.DataFrame({"inhib":ct["Inhibitory"],"other":ct.sum(axis=1)-ct["Inhibitory"]}).reindex(order)
fd=sst.fisher_exact(inh_tab.loc[[order[0],order[2]]].values); print("Fisher dominant vs rare inhibitory:",fd)
p_chi=sst.chi2_contingency(inh_tab.values)[1]; p_kw=sst.kruskal(*data).pvalue; print(p_chi,p_kw)
xl2=["Dominant\n≥1%","Interm.\n0.1–1%","Rare\n<0.1%"]
for a_ in axs: a_.set_xticklabels(xl2)
axs[1].set_title(f"Inhibitory share: 67% → 39%",loc="left")
axs[1].text(0.98,0.79,f"χ² p = {p_chi:.2f}",transform=axs[1].transAxes,ha="right",fontsize=6)
for t in list(axs[2].texts):
    if t.get_text()=="no effect": t.remove()
axs[2].text(0.55,4,"no effect",fontsize=6,color="#555555",ha="left",va="bottom")
axs[2].text(0.98,0.97,f"Kruskal–Wallis p = {p_kw:.2f}",transform=axs[2].transAxes,ha="right",va="top",fontsize=6)
axs[2].set_ylim(-215,150)
fig.savefig("AmphiBac_abundance_class_function.png",dpi=300,bbox_inches="tight")
host.view_image("AmphiBac_abundance_class_function.png")
summ=pd.DataFrame({"abundance_class":order,"total_sOTUs_in_class":[tot_cls[c] for c in order],"Bd_tested_database_taxa":n_t.values,
  "n_inhibitory":ct["Inhibitory"].values,"n_no_effect_mixed":ct["No effect / mixed"].values,"n_enhancing":ct["Enhancing"].values,
  "pct_inhibitory":pct["Inhibitory"].round(1).values,"pct_no_effect_mixed":pct["No effect / mixed"].round(1).values,"pct_enhancing":pct["Enhancing"].round(1).values,
  "taxa_with_score":[len(d) for d in data],"median_taxon_inhib_score":[round(np.median(d),1) for d in data],"IQR_low":[round(np.percentile(d,25),1) for d in data],"IQR_high":[round(np.percentile(d,75),1) for d in data],
  "isolates_scored":iso_cls.groupby("abundance_class").inhib_score.count().reindex(order).values,"median_isolate_inhib_score":iso_cls.groupby("abundance_class").inhib_score.median().reindex(order).round(1).values})
tests=pd.DataFrame([dict(test="Chi-square: inhibitory vs not, across 3 abundance classes",p=p_chi),dict(test="Fisher exact: inhibitory share, dominant vs rare",odds_ratio=fd.statistic,p=fd.pvalue),
  dict(test="Kruskal-Wallis: taxon median inhibition score across classes",p=p_kw)])
notes=pd.DataFrame({"note":["Taxa = amphibian meta-analysis sOTUs matched by >=1 Bd-tested AmphiBac 2025.1 isolate (exact 90-bp V4 containment, or <=3 mismatches for Sanger error).",
 "Abundance class from sOTU mean relative abundance across 2,411 meta-analysis samples (>=1,000 reads): dominant >=1%, intermediate 0.1-1%, rare <0.1%.",
 "Function per taxon: Inhibitory if >=50% of its Bd-tested isolates are inhibitory; Enhancing if >=50% facilitating; otherwise No effect / mixed.",
 "Bd inhibition score = (1 - Proportional.Growth.Bd) x 100, available for 96-well cell-free supernatant (CFS) assays only (agar-plate assays have no score). 100 = complete inhibition, 0 = no effect, negative = enhanced Bd growth. Panel c uses each taxon's median isolate score."]})
T_out=T.reset_index().sort_values("mean_rel_pct",ascending=False)
with pd.ExcelWriter("AmphiBac_abundance_class_function.xlsx") as w:
    summ.to_excel(w,sheet_name="summary_by_class",index=False); tests.to_excel(w,sheet_name="tests",index=False)
    T_out.to_excel(w,sheet_name="taxa",index=False); notes.to_excel(w,sheet_name="methods_notes",index=False)
print(summ.T.to_string())
