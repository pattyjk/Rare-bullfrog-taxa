# session cell 133 | 2026-10-18 11:37 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_008_AmphiBac_abundance_class_definitions.png","sha256":"e0c9b2f4a04dc42a321486747f3a58a5ec732591ebac39296a535d0c7da67993","preview_sha256":"e0c9b2f4a04dc42a321486747f3a58a5ec732591ebac39296a535d0c7da67993","preview_content_type":"image/png"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_abundance_class_definitions.png","sha256":"7b60ad895685866b553eca45e1b104763ec05f673794f284a501b972a844784f","preview_sha256":"f9b6c438daafadbbb3fbb245cc915a945b717eb8457c6edb28b0a0047326ea21","preview_content_type":"image/jpeg"}]


import matplotlib as mpl, matplotlib.pyplot as plt
apply_figure_style(sizes=(8,7,6))
defs=[("Global mean (all 2,411 samples, zeros included)",T.assign(cls=T.abundance_class),order,"taxa"),
      ("Mean when present (samples where detected)",T.assign(cls=T.class_present),order,"taxa"),
      ("Host-specific: mean on the isolate's source host",U.assign(cls=U.host_class),order4,"taxon × source-host units")]
labs={"Dominant (≥1%)":"Dom.\n≥1%","Intermediate (0.1–1%)":"Interm.\n0.1–1%","Rare (<0.1%)":"Rare\n<0.1%","Not detected on source host":"Not\ndetected"}
fcs=[("Inhibitory","#2166AC"),("No effect / mixed","#E08214"),("Enhancing","#7B3294")]
fig,axs=plt.subplots(3,3,figsize=(7.4,7.6),gridspec_kw={"wspace":0.5,"hspace":0.75,"width_ratios":[0.85,1,1.1]})
stats_rows=[]
for r,(title,D,ordr,unit) in enumerate(defs):
    n_=D.cls.value_counts().reindex(ordr).fillna(0).astype(int); x=np.arange(len(ordr))
    ax=axs[r,0]; ax.bar(x,n_.values,color="#4D4D4D",width=0.6)
    for i,v in enumerate(n_.values): ax.text(i,v+n_.max()*0.03,f"{v}",ha="center",va="bottom",fontsize=6)
    ax.set_ylim(0,n_.max()*1.2); ax.set_xticks(x); ax.set_xticklabels([labs[c] for c in ordr]); ax.set_ylabel(f"Bd-tested {unit.split(' ')[0]}\n(n)" if r<2 else "Taxon × host\nunits (n)")
    ct_=pd.crosstab(D.cls,D.function_class).reindex(ordr).fillna(0)[[f for f,_ in fcs]]; pc_=ct_.div(ct_.sum(axis=1),axis=0)*100
    ax2=axs[r,1]; b=np.zeros(len(ordr))
    for f,col in fcs:
        v=pc_[f].values; ax2.bar(x,v,bottom=b,color=col,width=0.6,label=f,edgecolor="white",lw=0.5)
        for i,vv in enumerate(v):
            if vv>=8: ax2.text(i,b[i]+vv/2,f"{vv:.0f}",ha="center",va="center",fontsize=6,color="white")
        b+=v
    det=ct_.loc[order]; p_inh=sst.chi2_contingency(np.c_[det.Inhibitory,det.sum(axis=1)-det.Inhibitory])[1]
    ax2.set_ylim(0,118); ax2.set_yticks([0,25,50,75,100]); ax2.set_xticks(x); ax2.set_xticklabels([labs[c] for c in ordr]); ax2.set_ylabel("Units by Bd\nfunction (%)")
    ax2.text(0.5,0.9,f"inhibitory vs not: χ² p = {p_inh:.2g}" if p_inh>=0.001 else "inhibitory vs not: χ² p < 0.001",transform=ax2.transAxes,ha="center",fontsize=6)
    dd=[D[(D.cls==c)&D.median_inhib_score.notna()].median_inhib_score.values for c in ordr]
    ax3=axs[r,2]; bp=ax3.boxplot(dd,widths=0.55,showfliers=False,patch_artist=True,medianprops=dict(color="black",lw=1.1))
    for p_ in bp["boxes"]: p_.set_facecolor("#D1E5F0"); p_.set_edgecolor("#333333")
    rng=np.random.default_rng(r)
    for i,d_ in enumerate(dd,1):
        ax3.scatter(i+rng.uniform(-0.17,0.17,len(d_)),np.clip(d_,-195,145),s=1.5,color="#2166AC",alpha=0.3,lw=0,zorder=3)
        ax3.text(i,-190,f"{len(d_)}",ha="center",fontsize=6)
    ax3.axhline(0,color="#666666",lw=0.6,ls="--"); ax3.set_ylim(-200,150)
    p_kw=sst.kruskal(*dd[:3]).pvalue
    ax3.text(0.98,0.98,f"KW p = {p_kw:.2g}" if p_kw>=0.001 else "KW p < 0.001",transform=ax3.transAxes,ha="right",va="top",fontsize=6)
    ax3.set_xticks(range(1,len(ordr)+1)); ax3.set_xticklabels([labs[c] for c in ordr]); ax3.set_ylabel("Inhibition score\n(median, %)")
    axs[r,0].set_title(title,loc="left",x=0,pad=14,fontsize=8,fontweight="bold") 
    medians=[round(np.median(d),1) if len(d) else np.nan for d in dd]
    for c,nv,pin,med in zip(ordr,n_.values,pc_.Inhibitory.values,medians): stats_rows.append(dict(definition=title,abundance_class=c,n_units=nv,pct_inhibitory=round(pin,1),median_score=med))
    stats_rows.append(dict(definition=title,abundance_class="TEST",chi2_inhibitory_p=p_inh,kruskal_p=p_kw))
axs[0,1].legend(frameon=False,loc="upper center",ncol=3,fontsize=6,handlelength=0.8,handletextpad=0.3,columnspacing=0.5,bbox_to_anchor=(0.5,1.32))
for i,a_ in enumerate(axs.flat): panel_letter(a_,"abcdefghi"[i])
fig.savefig("AmphiBac_abundance_class_definitions.png",dpi=300,bbox_inches="tight")
r_=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r_)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print([(a_.get_text()[:25],b_.get_text()[:25]) for i,(a_,ba) in enumerate(texts) for b_,bb in texts[i+1:] if ba.overlaps(bb)][:12])
host.view_image("AmphiBac_abundance_class_definitions.png")
