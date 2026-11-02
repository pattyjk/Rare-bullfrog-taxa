# session cell 121 | 2026-10-18 11:30 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_001_AmphiBac_abundance_class_function.png","sha256":"9dbe2fcfdfdefb62d79f1847fda37ea7f8cc15e396eaf1edde690184079ad5d5","preview_sha256":"9dbe2fcfdfdefb62d79f1847fda37ea7f8cc15e396eaf1edde690184079ad5d5","preview_content_type":"image/png"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_abundance_class_function.png","sha256":"bf7215399016bfc21ab229b3ea47afda3e495a01c9be98c11fd7ef9aae6a9aee","preview_sha256":"c9d454e5a81c5de433bac8d4b8009a84f43ee7defae5de570f072975bc56a9e6","preview_content_type":"image/jpeg"}]


import matplotlib as mpl, matplotlib.pyplot as plt
apply_figure_style(sizes=(8,7,6))
xl=["Dominant\n≥1%","Intermediate\n0.1–1%","Rare\n<0.1%"]
tot_cls={"Dominant (≥1%)":15,"Intermediate (0.1–1%)":82,"Rare (<0.1%)":48334}
n_t=T.abundance_class.value_counts().reindex(order)
fig,axs=plt.subplots(1,3,figsize=(7.4,2.9),gridspec_kw={"wspace":0.5,"width_ratios":[0.85,0.95,1.1]})
ax=axs[0]; ax.bar(range(3),n_t.values,color="#4D4D4D",width=0.6)
for i,c in enumerate(order):
    lab=f"{n_t[c]}\nof {tot_cls[c]/1000:.1f}k" if tot_cls[c]>=1000 else f"{n_t[c]}\nof {tot_cls[c]}"
    ax.text(i,n_t[c]+8,lab,ha="center",va="bottom",fontsize=6)
ax.set_xticks(range(3)); ax.set_xticklabels(xl); ax.set_ylim(0,560); ax.set_ylabel("Bd-tested database taxa (sOTUs)")
ax.set_title("Most tested taxa are rare",loc="left")
ax2=axs[1]; fcs=[("Inhibitory","#2166AC"),("No effect / mixed","#E08214"),("Enhancing","#7B3294")]
pct=ct.div(ct.sum(1),axis=0)*100; bottom=np.zeros(3)
for fc,col in fcs:
    v=pct[fc].values; ax2.bar(range(3),v,bottom=bottom,color=col,width=0.6,label=fc,edgecolor="white",lw=0.5)
    for i,(vv,b,n_) in enumerate(zip(v,bottom,ct[fc].values)):
        if vv>=6: ax2.text(i,b+vv/2,f"{vv:.0f}%",ha="center",va="center",fontsize=6,color="white")
    bottom+=v
ax2.set_xticks(range(3)); ax2.set_xticklabels(xl); ax2.set_ylim(0,125); ax2.set_yticks([0,25,50,75,100])
ax2.set_ylabel("Taxa by Bd function (%)"); ax2.set_title("Share inhibitory is similar",loc="left")
ax2.legend(frameon=False,loc="upper center",ncol=3,fontsize=6,handlelength=0.9,columnspacing=0.6,borderaxespad=0,bbox_to_anchor=(0.5,1.0))
ax3=axs[2]
data=[T[(T.abundance_class==c)&T.median_inhib_score.notna()].median_inhib_score.values for c in order]
bp=ax3.boxplot(data,widths=0.55,showfliers=False,patch_artist=True,medianprops=dict(color="black",lw=1.2))
for p_ in bp["boxes"]: p_.set_facecolor("#D1E5F0"); p_.set_edgecolor("#333333")
rng=np.random.default_rng(2)
for i,d_ in enumerate(data,1): ax3.scatter(i+rng.uniform(-0.17,0.17,len(d_)),d_,s=3,color="#2166AC",alpha=0.45,lw=0,zorder=3)
ax3.axhline(0,color="#666666",lw=0.6,ls="--"); ax3.text(3.45,3,"no effect",fontsize=6,color="#555555",ha="right",va="bottom")
for i,d_ in enumerate(data,1): ax3.text(i,-205,f"n={len(d_)}",ha="center",fontsize=6)
ax3.set_xticks([1,2,3]); ax3.set_xticklabels(xl); ax3.set_ylim(-215,130)
ax3.set_ylabel("Bd inhibition score\n(taxon median, %)"); ax3.set_title("Inhibition strength is similar",loc="left")
for a_,l_ in zip(axs,"abc"): panel_letter(a_,l_)
fig.savefig("AmphiBac_abundance_class_function.png",dpi=300,bbox_inches="tight")
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print([(a_.get_text()[:25],b_.get_text()[:25]) for i,(a_,ba) in enumerate(texts) for b_,bb in texts[i+1:] if ba.overlaps(bb)][:10])
host.view_image("AmphiBac_abundance_class_function.png")
