# session cell 104 | 2026-10-17 18:13 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_cultured_vs_abundance.png","sha256":"5a8bb7a31ab4c81b73b40fe2ac4215092d279ad56df97278e7652fd2c9ef29ab","preview_sha256":"d055f86ffbb6afd22ab3686f2c70fa5e73fb326603028c5a34cc91f7f78d98eb","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_007_AmphiBac_cultured_vs_abundance.png","sha256":"d54b643f6ecc766531caeca1ed7ae1c9a66d1d5b912a6f8d0daab7a3fab91c14","preview_sha256":"d54b643f6ecc766531caeca1ed7ae1c9a66d1d5b912a6f8d0daab7a3fab91c14","preview_content_type":"image/png"}]


apply_figure_style(sizes=(8,7,6))
fig,axs=plt.subplots(1,3,figsize=(7.4,2.9),gridspec_kw={"wspace":0.55,"width_ratios":[1.15,0.85,1.1]})
ax=axs[0]; c_=cf.loc[order]; x=np.arange(3); w=0.38
ax.bar(x-w/2,c_.pct_cultured_exact,w,color="#2166AC",label="Exact V4 match")
ax.bar(x+w/2,c_.pct_cultured_3mm,w,color="#92C5DE",label="≤3 mismatches (≥96.7%)")
for i,r_ in enumerate(c_.itertuples()):
    if i<2:
        ax.text(i-w/2,r_.pct_cultured_exact+2,f"{int(r_.n_cultured_exact)}/{int(r_.n_sOTUs)}",ha="center",va="bottom",fontsize=6)
        ax.text(i+w/2,r_.pct_cultured_3mm+2,f"{int(r_.n_cultured_3mm)}",ha="center",va="bottom",fontsize=6)
ax.text(2,4,f"{int(c_.iloc[2].n_cultured_exact)}/48.3k\n{int(c_.iloc[2].n_cultured_3mm)}",ha="center",va="bottom",fontsize=6)
ax.set_xticks(x); ax.set_xticklabels(["Dominant\n≥1%","Interm.\n0.1–1%","Rare\n<0.1%"])
ax.set_ylim(0,125); ax.set_yticks([0,20,40,60,80,100]); ax.set_ylabel("sOTUs with a cultured\nAmphiBac match (%)")
ax.set_title("All dominant sOTUs are cultured",loc="left"); ax.legend(frameon=False,loc="upper right",fontsize=6,borderaxespad=0)
ax2=axs[1]; d1_=PSm.pct_sOTUs_cultured_3mm.values; d2_=PSm.pct_reads_cultured_3mm.values
bp=ax2.boxplot([d1_,d2_],widths=0.5,showfliers=False,patch_artist=True,medianprops=dict(color="black",lw=1.2))
for p_,col in zip(bp["boxes"],["#D9D9D9","#92C5DE"]): p_.set_facecolor(col); p_.set_edgecolor("#333333")
rng=np.random.default_rng(1)
for i,d_ in enumerate([d1_,d2_],1):
    ax2.scatter(i+rng.uniform(-0.18,0.18,len(d_)),d_,s=1,color="#333333",alpha=0.12,lw=0,zorder=3)
    ax2.text(i+0.29,np.median(d_),f"{np.median(d_):.0f}%",va="center",fontsize=6)
ax2.set_xticks([1,2]); ax2.set_xticklabels(["sOTUs","Reads"]); ax2.set_ylim(-2,102)
ax2.set_ylabel("Share matched to cultured\nisolates per sample (%)"); ax2.set_title("Few sOTUs, most reads",loc="left")
ax3=axs[2]; xm=PSm.cultured_sOTUs_3mm.values; yr=PSm.pct_reads_cultured_3mm.values
ax3.scatter(xm,yr,s=2,color="#2166AC",alpha=0.25,lw=0)
bins=np.arange(0,xm.max()+10,10); bc=[];bm=[]
for lo in bins[:-1]:
    sel=(xm>=lo)&(xm<lo+10)
    if sel.sum()>=20: bc.append(lo+5); bm.append(np.median(yr[sel]))
ax3.plot(bc,bm,color="black",lw=1.2,marker="o",ms=2.5,label="Median per 10-match bin")
rho=sst.spearmanr(xm,yr).statistic
ax3.text(0.97,0.04,f"Spearman ρ = {rho:.2f}\nn = {len(xm):,} samples",transform=ax3.transAxes,ha="right",va="bottom",fontsize=6)
ax3.set_xlabel("Cultured sOTUs per sample (n)"); ax3.set_ylabel("Reads matched to cultured\nisolates (%)")
ax3.set_ylim(-2,102); ax3.margins(x=0.03)
ax3.set_title("Even a few matches carry most reads",loc="left"); ax3.legend(frameon=False,loc="upper right",fontsize=6,borderaxespad=0)
for a_,l_ in zip(axs,"abc"): panel_letter(a_,l_)
fig.savefig("AmphiBac_cultured_vs_abundance.png",dpi=300,bbox_inches="tight")
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print([(a_.get_text(),b_.get_text()) for i,(a_,ba) in enumerate(texts) for b_,bb in texts[i+1:] if ba.overlaps(bb)][:8])
print("median %reads for samples with <=10 matches:",round(np.median(yr[xm<=10]),1),"n",int((xm<=10).sum()))
host.view_image("AmphiBac_cultured_vs_abundance.png")
