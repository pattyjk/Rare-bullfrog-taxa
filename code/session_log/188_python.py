# session cell 188 | 2026-11-02 11:15 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_multihit_sensitivity.png","sha256":"227eedc1b6b9c6b33cda85c50b85bd88c26c1884ef6471b0f16cb392d0458114","preview_sha256":"c9942feadffb66cdbfd16e51b8e533a6eeb8e4995bcb9208e285b4b79809a70f","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_001_AmphiBac_multihit_sensitivity.png","sha256":"a1ea2ed98850731ef064d6e31c770ec541279ec1fa55335a1cee3ff26959f819","preview_sha256":"a1ea2ed98850731ef064d6e31c770ec541279ec1fa55335a1cee3ff26959f819","preview_content_type":"image/png"}]

import matplotlib as mpl, matplotlib.pyplot as plt
apply_figure_style(sizes=(8,7,6))
colc={'Dominant':'#08306B','Intermediate':'#4292C6','Rare':'#A6A6A6'}
short={'Tied nearest -> most abundant sOTU (used so far)':'Most abundant tied sOTU (used so far)',
 'Tied nearest -> rarest sOTU':'Rarest tied sOTU','Tied nearest -> all tied sOTUs':'All tied sOTUs',
 'Drop isolates with tied nearest sOTUs':'Drop tied isolates','Exact matches only (0 mismatches)':'Exact matches only',
 'Tied nearest -> random sOTU (median of 500 draws)':'Random tied sOTU (500 draws)',
 'Majority of isolates inhibitory (used so far)':'Majority inhibitory (used so far)',
 'Any isolate inhibitory (AmphiBac subset-mode logic)':'Any isolate inhibitory',
 'All isolates inhibitory':'All isolates inhibitory',
 'One random isolate per sOTU (top-hit analogue; median of 500 draws)':'One random isolate (top-hit; 500 draws)',
 'Unanimous sOTUs only (conflicting sOTUs dropped)':'Unanimous sOTUs only',
 'Majority, sOTUs with >=2 isolates only':'Majority, sOTUs with \u22652 isolates'}
fig,axs=plt.subplots(2,1,figsize=(6.6,4.6),sharex=True,gridspec_kw={'height_ratios':[6,6],'hspace':0.35})
for ax,grp,rand,title in [(axs[0],'A. isolate-to-sOTU tie rule (majority label)',rr,'Isolates tied between sOTUs: how they are assigned'),
                          (axs[1],'B. sOTU label rule (original tie rule)',rr2,'sOTUs with conflicting isolate labels: how the sOTU is called')]:
    d=R_[R_.group==grp].reset_index(drop=True)
    yy=np.arange(len(d))[::-1]
    for y,(_,r) in zip(yy,d.iterrows()):
        ax.plot([min(r.inh_D,r.inh_I,r.inh_R),max(r.inh_D,r.inh_I,r.inh_R)],[y,y],color='#DDDDDD',lw=1,zorder=1)
        israndom='500 draws' in r.scenario
        for c,k in [('Rare','inh_R'),('Intermediate','inh_I'),('Dominant','inh_D')]:
            if israndom:
                lo,hi=rand[k].quantile([.025,.975]); ax.plot([lo,hi],[y,y],color=colc[c],lw=2.2,alpha=0.45,solid_capstyle='butt',zorder=2)
            ax.scatter(r[k],y,s=26,color=colc[c],zorder=3,edgecolor='white',lw=0.4)
        p=r.chi2_p; ax.text(101,y,(f"p = {p:.2g}" if p>=0.001 else "p < 0.001")+(" (median)" if israndom else ""),va='center',fontsize=6,color='#333333')
    ax.set_yticks(yy); ax.set_yticklabels([short[s] for s in d.scenario])
    for t in ax.get_yticklabels():
        if 'used so far' in t.get_text(): t.set_fontweight('bold')
    ax.set_ylim(-0.6,len(d)-0.4); ax.set_title(title,loc='left')
axs[1].set_xlim(0,100); axs[1].set_xlabel('sOTUs called Bd-inhibitory (%)')
from matplotlib.lines import Line2D
h=[Line2D([],[],marker='o',ls='',color=colc[c],label=l,markersize=5) for c,l in [('Dominant','Dominant \u22651%'),('Intermediate','Intermediate 0.1\u20131%'),('Rare','Rare <0.1%')]]
axs[0].legend(handles=h,frameon=False,loc='lower left',bbox_to_anchor=(0.0,1.12),ncol=3,fontsize=6,handletextpad=0.2,columnspacing=0.8,borderaxespad=0)
axs[0].text(1.0,1.16,'\u03c7\u00b2 test, class \u00d7 inhibitory',transform=axs[0].transAxes,ha='right',fontsize=6,color='#333333')
for a_,l_ in zip(axs,'ab'): panel_letter(a_,l_)
fig.savefig('AmphiBac_multihit_sensitivity.png',dpi=300,bbox_inches='tight')
r_=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r_)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print([(a.get_text()[:25],b.get_text()[:25]) for i,(a,ba) in enumerate(texts) for b,bb in texts[i+1:] if ba.overlaps(bb)])
host.view_image('AmphiBac_multihit_sensitivity.png')