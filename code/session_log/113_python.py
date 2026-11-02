# session cell 113 | 2026-10-17 21:34 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_source_host_abundance.png","sha256":"3dd56ccf0569ca5d86203d7407d6ae12fb0baf15dc0d0810b68ab713a64fc078","preview_sha256":"26957528ae41ec215bc89cc6ca3ce02f60b12ae919c8d2fa8f5e0a377d33bbce","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_001_AmphiBac_source_host_abundance.png","sha256":"18cfb1232aa209dc2ff76dd4551e3499c65a6bff2132d4ca9e43bb1601592daa","preview_sha256":"18cfb1232aa209dc2ff76dd4551e3499c65a6bff2132d4ca9e43bb1601592daa","preview_content_type":"image/png"}]


import matplotlib as mpl, matplotlib.pyplot as plt
apply_figure_style(sizes=(8,7,6))
S=R2[R2.home_group_mean.notna()].copy()
f=lambda v:np.log10(v+pc)
fig,axs=plt.subplots(1,2,figsize=(7.0,3.2),gridspec_kw={"wspace":0.4,"width_ratios":[1.1,1]})
ax=axs[0]; cols={"same host species":"#2166AC","same host genus":"#67A9CF","same host family":"#BDBDBD"}
for lv in ["same host family","same host genus","same host species"]:
    s=S[S.home_level==lv]; ax.scatter(f(s.home_group_mean),f(s.median_other_groups),s=4,alpha=0.45,lw=0,color=cols[lv],label=f"{lv.replace('same host ','Source host ').capitalize() if False else lv.capitalize()} (n={len(s):,})")
lim=(f(0)-0.15,1.3); ax.plot(lim,lim,color="black",lw=0.7,ls="--"); ax.set_xlim(lim); ax.set_ylim(lim)
tk=[0,1e-3,1e-2,1e-1,1,10]; tl=["0","0.001","0.01","0.1","1","10"]
ax.set_xticks([f(t) for t in tk]); ax.set_xticklabels(tl); ax.set_yticks([f(t) for t in tk]); ax.set_yticklabels(tl)
ax.set_xlabel("Mean relative abundance on source host (%)"); ax.set_ylabel("Median across other host\nspecies/genera/families (%)")
ax.text(0.97,0.05,"above line: higher on other hosts",transform=ax.transAxes,ha="right",fontsize=6,color="#444444")
ax.set_title("Cultured taxa are only slightly enriched on their source host",loc="left",fontsize=8)
ax.legend(frameon=False,loc="upper left",fontsize=6,markerscale=2.5,handletextpad=0.2,borderaxespad=0)
ax2=axs[1]; sp=S[S.home_level=="same host species"]
bins_=np.linspace(0,100,11)
for fc,col in [("Inhibitory","#2166AC"),("Non-inhibitory/facilitating","#E08214")]:
    v=sp[sp.function==fc].pct_other_groups_below_home
    ax2.hist(v,bins=bins_,histtype="step",lw=1.4,color=col,density=False,label=f"{fc} (n={len(v)}, median {v.median():.0f}%)")
    ax2.axvline(v.median(),color=col,lw=0.8,ls=":")
ax2.set_xlabel("Other host species with lower abundance\nthan the source host (%)"); ax2.set_ylabel("Isolate–host units (n)")
ax2.set_title("Source host ranks mid-table, for both functions",loc="left",fontsize=8)
ax2.legend(frameon=False,loc="upper center",fontsize=6,bbox_to_anchor=(0.5,1.0))
ax2.set_ylim(0,ax2.get_ylim()[1]*1.25)
panel_letter(axs[0],"a"); panel_letter(axs[1],"b")
fig.savefig("AmphiBac_source_host_abundance.png",dpi=300,bbox_inches="tight")
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print([(a_.get_text()[:30],b_.get_text()[:30]) for i,(a_,ba) in enumerate(texts) for b_,bb in texts[i+1:] if ba.overlaps(bb)][:8])
print(sp.groupby("function").pct_other_groups_below_home.median().to_dict())
host.view_image("AmphiBac_source_host_abundance.png")
