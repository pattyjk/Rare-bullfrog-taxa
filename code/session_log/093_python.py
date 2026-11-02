# session cell 93 | 2026-10-17 17:53 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_rank_abundance.png","sha256":"d207a5b2cadfbbcbfe4cd085f6720c00edd5b038a00bf8cbaf875a67a12febb9","preview_sha256":"8c92e1c60ac6d8d77721dce26aaf62aa4276d0d5419b77b2631ea1defc1d0df5","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_001_AmphiBac_rank_abundance.png","sha256":"8095d121cf4ef49ec3bffb704eeba5189e227c2e1fb3b65c37617896ea5f691c","preview_sha256":"8095d121cf4ef49ec3bffb704eeba5189e227c2e1fb3b65c37617896ea5f691c","preview_content_type":"image/png"}]


import matplotlib as mpl, matplotlib.pyplot as plt
apply_figure_style(sizes=(8,7,6))
P=present.sort_values("mean_rel_pct",ascending=False).copy(); P["rank"]=np.arange(1,len(P)+1)
fig,axs=plt.subplots(1,2,figsize=(7.2,3.3),gridspec_kw={"width_ratios":[2.1,1],"wspace":0.35})
ax=axs[0]
ax.plot(P["rank"],P.mean_rel_pct,color="#BBBBBB",lw=1.2,zorder=1,label=f"All amphibian sOTUs (n={len(P)/1000:.1f}k)")
colf={"Inhibitory":"#2166AC","Non-inhibitory/mixed":"#E08214","Facilitating":"#7B3294"}
for fc,col in colf.items():
    s=agg[agg.function_class==fc]; ax.scatter(P.loc[s.index,"rank"],P.loc[s.index,"mean_rel_pct"],s=9,color=col,alpha=0.75,lw=0,zorder=2,label=f"AmphiBac {fc.lower()} (n={len(s)})")
ax.axhline(0.1,color="#666666",lw=0.6,ls=":",zorder=0); ax.text(P["rank"].max()*0.9,0.115,"0.1% (rare below)",ha="right",va="bottom",fontsize=6,color="#555555")
labs={"BF3":"BF3","BF6":"BF6","BF5":"BF5","BF1":"BF1"}
offs={"BF3":(-1.0,-1.5),"BF6":(0.15,1.0),"BF5":(-0.2,0.6),"BF1":(0.5,0.4)}
for _,r_ in I2[(I2.note=="")&(I2.meta_mean_rel_pct>0)].iterrows():
    s90=r_.V4_120bp[:90]; x,y=P.loc[s90,"rank"],P.loc[s90,"mean_rel_pct"]
    ax.scatter([x],[y],marker="*",s=70,color="black",zorder=4,edgecolor="white",lw=0.4)
    dx,dy=offs[r_.isolate]; tx,ty=x*10**dx,y*10**dy
    ax.annotate(r_.isolate,(x,y),(tx,ty),fontsize=7,arrowprops=dict(arrowstyle="-",lw=0.5,color="black",shrinkA=0,shrinkB=3),ha="center",va="center")
ax.scatter([],[],marker="*",s=70,color="black",label="Your isolates (genome V4)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("sOTU abundance rank"); ax.set_ylabel("Mean relative abundance (%)")
ax.set_title("Cultured isolates match the most abundant skin sOTUs",loc="left")
ax.legend(frameon=False,loc="lower left",fontsize=6,handletextpad=0.3,borderaxespad=0.2)
ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v,_:f"{v:g}"))
ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v,_:{1:"1",10:"10",100:"100",1000:"1k",10000:"10k",100000:"100k"}.get(round(v),"")))
ax.margins(0.04)
ax2=axs[1]
cl=cls_pct.set_index("abundance_class").loc[["dominant (>=1%)","intermediate (0.1-1%)","rare (<0.1%)"]]
xl=["Dominant\n≥1%","Intermediate\n0.1–1%","Rare\n<0.1%"]
ax2.bar(range(3),cl.pct_inhibitory,color="#2166AC",width=0.65)
for i,(p_,n_,k_) in enumerate(zip(cl.pct_inhibitory,cl.n_sOTUs,cl.n_inhibitory)): ax2.text(i,p_+1.5,f"{int(k_)}/{int(n_)}",ha="center",va="bottom",fontsize=6)
ax2.set_xticks(range(3)); ax2.set_xticklabels(xl); ax2.set_ylim(0,72)
ax2.set_ylabel("Bd-inhibitory sOTUs (%)"); ax2.set_title("Rare sOTUs are not more often inhibitory",loc="left")
panel_letter(axs[0],"a"); panel_letter(axs[1],"b")
fig.savefig("AmphiBac_rank_abundance.png",dpi=300,bbox_inches="tight")
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print([(a_.get_text(),b_.get_text()) for i,(a_,ba) in enumerate(texts) for b_,bb in texts[i+1:] if ba.overlaps(bb)][:8])
host.view_image("AmphiBac_rank_abundance.png")
