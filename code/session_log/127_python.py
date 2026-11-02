# session cell 127 | 2026-10-18 11:33 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_006_AmphiBac_abundance_class_function_isolates.png","sha256":"373004bb5ee81d454b4b1a835a0ab8d065f06b1e0a8b55255f967a8d5d2ab63f","preview_sha256":"373004bb5ee81d454b4b1a835a0ab8d065f06b1e0a8b55255f967a8d5d2ab63f","preview_content_type":"image/png"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_abundance_class_function_isolates.png","sha256":"223cd0809a61d60f14e1baabf3aa781a036f3ce9fe22c4d7f61818347ab8951f","preview_sha256":"3736b83798f9e21a01cca20a1f92c074e4498ac759c4b46b976dc29b223ec688","preview_content_type":"image/jpeg"}]


allsc=np.concatenate(di); lo_,hi_=-215,150; n_out=int(((allsc<lo_)|(allsc>hi_)).sum()); print("outside view:",n_out, allsc.min().round(0), allsc.max().round(0))
apply_figure_style(sizes=(8,7,6))
xl2=["Dominant\n≥1%","Interm.\n0.1–1%","Rare\n<0.1%"]
fig,axs=plt.subplots(1,3,figsize=(7.4,2.9),gridspec_kw={"wspace":0.5,"width_ratios":[0.85,0.95,1.1]})
ax=axs[0]; ax.bar(range(3),n_i.values,color="#4D4D4D",width=0.6)
for i,c in enumerate(order): ax.text(i,n_i[c]+60,f"{n_i[c]:,}",ha="center",va="bottom",fontsize=6)
ax.set_xticks(range(3)); ax.set_xticklabels(xl2); ax.set_ylim(0,3700); ax.set_ylabel("Bd-tested database isolates (n)")
ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v,_:f"{v/1000:g}k" if v else "0"))
ax.set_title("Over half of isolates are rare taxa",loc="left")
ax2=axs[1]; bottom=np.zeros(3)
for fc,col in [("Inhibitory","#2166AC"),("No effect","#E08214"),("Enhancing","#7B3294")]:
    v=pcti[fc].values; ax2.bar(range(3),v,bottom=bottom,color=col,width=0.6,label=fc,edgecolor="white",lw=0.5)
    for i,vv in enumerate(v):
        if vv>=6: ax2.text(i,bottom[i]+vv/2,f"{vv:.0f}%",ha="center",va="center",fontsize=6,color="white")
    bottom+=v
ax2.set_xticks(range(3)); ax2.set_xticklabels(xl2); ax2.set_ylim(0,125); ax2.set_yticks([0,25,50,75,100])
ax2.set_ylabel("Isolates by Bd function (%)"); ax2.set_title("Inhibitory share: 68% → 40%",loc="left")
ax2.legend(frameon=False,loc="upper center",ncol=3,fontsize=6,handlelength=0.8,handletextpad=0.3,columnspacing=0.5,borderaxespad=0,bbox_to_anchor=(0.52,1.0))
ax2.text(0.5,0.855,"χ² p < 0.001",transform=ax2.transAxes,ha="center",fontsize=6)
ax3=axs[2]
bp=ax3.boxplot(di,widths=0.55,showfliers=False,patch_artist=True,medianprops=dict(color="black",lw=1.2))
for p_ in bp["boxes"]: p_.set_facecolor("#D1E5F0"); p_.set_edgecolor("#333333")
rng=np.random.default_rng(3)
for i,d_ in enumerate(di,1):
    dd=np.clip(d_,lo_+5,hi_-5); ax3.scatter(i+rng.uniform(-0.18,0.18,len(dd)),dd,s=1.2,color="#2166AC",alpha=0.18,lw=0,zorder=3)
    ax3.text(i,-205,f"n={len(d_):,}",ha="center",fontsize=6)
ax3.axhline(0,color="#666666",lw=0.6,ls="--"); ax3.text(0.55,4,"no effect",fontsize=6,color="#555555",va="bottom")
ax3.text(0.98,0.97,"Kruskal–Wallis p < 0.001",transform=ax3.transAxes,ha="right",va="top",fontsize=6)
ax3.set_xticks([1,2,3]); ax3.set_xticklabels(xl2); ax3.set_ylim(lo_,hi_)
ax3.set_ylabel("Bd inhibition score\n(per isolate, %)"); ax3.set_title("Inhibition is stronger in dominant taxa",loc="left")
for a_,l_ in zip(axs,"abc"): panel_letter(a_,l_)
fig.savefig("AmphiBac_abundance_class_function_isolates.png",dpi=300,bbox_inches="tight")
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print([(a_.get_text()[:25],b_.get_text()[:25]) for i,(a_,ba) in enumerate(texts) for b_,bb in texts[i+1:] if ba.overlaps(bb)][:10])
host.view_image("AmphiBac_abundance_class_function_isolates.png")
