# session cell 98 | 2026-10-17 17:59 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_cultured_vs_abundance.png","sha256":"e2bcc26197260985a2e9d59b4a606485dd256b5258de99b6b8d67037243e7993","preview_sha256":"d6c9c74599391eb455886da1577797dbcb2a89bfe5be169153bcb08b24fdfcc7","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_003_AmphiBac_cultured_vs_abundance.png","sha256":"6940989de59f7b5eaa12e7f6ce8ea76cf83368c86b76a4f87d1235e41d985362","preview_sha256":"6940989de59f7b5eaa12e7f6ce8ea76cf83368c86b76a4f87d1235e41d985362","preview_content_type":"image/png"}]


apply_figure_style(sizes=(8,7,6))
fig,axs=plt.subplots(1,2,figsize=(6.8,3.0),gridspec_kw={"wspace":0.45})
ax=axs[0]; order=["dominant (>=1%)","intermediate (0.1-1%)","rare (<0.1%)"]; c_=cf.loc[order]
x=np.arange(3); w=0.38
ax.bar(x-w/2,c_.pct_cultured_exact,w,color="#2166AC",label="Exact V4 match")
ax.bar(x+w/2,c_.pct_cultured_99,w,color="#92C5DE",label="≤1 mismatch (≥98.9%)")
for i,r_ in enumerate(c_.itertuples()):
    ax.text(i-w/2,r_.pct_cultured_exact+2,f"{int(r_.n_cultured_exact)}/{int(r_.n_sOTUs):,}",ha="center",va="bottom",fontsize=6)
    ax.text(i+w/2,r_.pct_cultured_99+2,f"{int(r_.n_cultured_99)}",ha="center",va="bottom",fontsize=6)
ax.set_xticks(x); ax.set_xticklabels(["Dominant\n≥1%","Intermediate\n0.1–1%","Rare\n<0.1%"])
ax.set_ylim(0,112); ax.set_ylabel("sOTUs with a cultured\nAmphiBac match (%)")
ax.set_title("Every dominant sOTU has been cultured; <1% of rare ones",loc="left")
ax.legend(frameon=False,loc="upper right",fontsize=6)
ax2=axs[1]
d1_=PSm.pct_sOTUs_cultured_exact.values; d2_=PSm.pct_reads_cultured_exact.values
bp=ax2.boxplot([d1_,d2_],widths=0.5,showfliers=False,patch_artist=True,medianprops=dict(color="black",lw=1.2))
for p_,col in zip(bp["boxes"],["#D9D9D9","#2166AC"]): p_.set_facecolor(col); p_.set_edgecolor("#333333")
rng=np.random.default_rng(1)
for i,d_ in enumerate([d1_,d2_],1): ax2.scatter(i+rng.uniform(-0.18,0.18,len(d_)),d_,s=1,color="#333333",alpha=0.12,lw=0,zorder=3)
for i,d_ in enumerate([d1_,d2_],1): ax2.text(i+0.3,np.median(d_),f"{np.median(d_):.0f}%",va="center",fontsize=6)
ax2.set_xticks([1,2]); ax2.set_xticklabels(["sOTUs\n(richness)","Reads"]); ax2.set_ylim(-2,102)
ax2.set_ylabel("Share matched to\ncultured isolates (%)")
ax2.set_title(f"Cultured taxa: few sOTUs, most reads",loc="left")
panel_letter(axs[0],"a"); panel_letter(axs[1],"b")
fig.savefig("AmphiBac_cultured_vs_abundance.png",dpi=300,bbox_inches="tight")
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print([(a_.get_text(),b_.get_text()) for i,(a_,ba) in enumerate(texts) for b_,bb in texts[i+1:] if ba.overlaps(bb)][:8], len(PSm))
host.view_image("AmphiBac_cultured_vs_abundance.png")
