# session cell 102 | 2026-10-17 18:02 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_006_AmphiBac_cultured_vs_abundance.png","sha256":"980ef41208765ba368db5aa3913905d4e2415e6894e921d826eaf48e53e28bcf","preview_sha256":"980ef41208765ba368db5aa3913905d4e2415e6894e921d826eaf48e53e28bcf","preview_content_type":"image/png"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_cultured_vs_abundance.png","sha256":"95ed858a756fd43bbb51a937998b11064e9a4e40d1b7d713336adcb2cdb5c3ed","preview_sha256":"0edec63d151f2859287334e3a5ae99e48bbb0ca4fe06e21237e5946200a90be4","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_near_matches.xlsx","sha256":"8dadb6841b24d92445f49b6af3586dc459e48bc3bc6ed8232ce72879e694354c"}]


apply_figure_style(sizes=(8,7,6))
fig,axs=plt.subplots(1,2,figsize=(6.8,3.0),gridspec_kw={"wspace":0.45})
ax=axs[0]; c_=cf.loc[order]; x=np.arange(3); w=0.38
ax.bar(x-w/2,c_.pct_cultured_exact,w,color="#2166AC",label="Exact V4 match")
ax.bar(x+w/2,c_.pct_cultured_3mm,w,color="#92C5DE",label="≤3 mismatches (≥96.7%)")
for i,r_ in enumerate(c_.itertuples()):
    if i<2:
        ax.text(i-w/2,r_.pct_cultured_exact+2,f"{int(r_.n_cultured_exact)}/{int(r_.n_sOTUs)}",ha="center",va="bottom",fontsize=6)
        ax.text(i+w/2,r_.pct_cultured_3mm+2,f"{int(r_.n_cultured_3mm)}",ha="center",va="bottom",fontsize=6)
ax.text(2,4,f"{int(c_.iloc[2].n_cultured_exact)}/48.3k ({c_.iloc[2].pct_cultured_exact:.1f}%)\n{int(c_.iloc[2].n_cultured_3mm)} ({c_.iloc[2].pct_cultured_3mm:.1f}%)",ha="center",va="bottom",fontsize=6)
ax.set_xticks(x); ax.set_xticklabels(["Dominant\n≥1%","Intermediate\n0.1–1%","Rare\n<0.1%"])
ax.set_ylim(0,112); ax.set_ylabel("sOTUs with a cultured\nAmphiBac match (%)")
ax.set_title("All dominant sOTUs have been cultured",loc="left"); ax.legend(frameon=False,loc="upper right",fontsize=6)
ax2=axs[1]
d1_=PSm.pct_sOTUs_cultured_3mm.values; d2_=PSm.pct_reads_cultured_3mm.values
bp=ax2.boxplot([d1_,d2_],widths=0.5,showfliers=False,patch_artist=True,medianprops=dict(color="black",lw=1.2))
for p_,col in zip(bp["boxes"],["#D9D9D9","#92C5DE"]): p_.set_facecolor(col); p_.set_edgecolor("#333333")
rng=np.random.default_rng(1)
for i,d_ in enumerate([d1_,d2_],1):
    ax2.scatter(i+rng.uniform(-0.18,0.18,len(d_)),d_,s=1,color="#333333",alpha=0.12,lw=0,zorder=3)
    ax2.text(i+0.3,np.median(d_),f"{np.median(d_):.0f}%",va="center",fontsize=6)
ax2.set_xticks([1,2]); ax2.set_xticklabels(["sOTUs\n(richness)","Reads"]); ax2.set_ylim(-2,102)
ax2.set_ylabel("Share matched to cultured\nisolates, ≤3 mismatches (%)"); ax2.set_title("Few sOTUs, most reads",loc="left")
panel_letter(axs[0],"a"); panel_letter(axs[1],"b")
fig.savefig("AmphiBac_cultured_vs_abundance.png",dpi=300,bbox_inches="tight")
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print([(a_.get_text(),b_.get_text()) for i,(a_,ba) in enumerate(texts) for b_,bb in texts[i+1:] if ba.overlaps(bb)][:8])
host.view_image("AmphiBac_cultured_vs_abundance.png")
mmsum["isolates_still_undetected_with_Bd_label"]=[int((~kn1.FastaID.isin(set(bt.FastaID[bt.detection.str.startswith('detected')])|set(nd.FastaID[nd.min_mismatch<=t]))).sum()) if t>0 else int((kn1.detection!='detected').sum()) for t in [0,1,2,3,5]]
with pd.ExcelWriter("AmphiBac_near_matches.xlsx") as w_:
    mmsum.to_excel(w_,sheet_name="threshold_summary",index=False); nd_out.to_excel(w_,sheet_name="unmatched_V4_isolates",index=False)
    cf.reset_index().rename(columns={"cls":"abundance_class"}).to_excel(w_,sheet_name="cultured_by_class",index=False)
    PSm.reset_index().rename(columns={"index":"SampleID"}).to_excel(w_,sheet_name="per_sample_cultured_share",index=False)
