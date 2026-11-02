# session cell 71 | 2026-10-15 22:05 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/antiSMASH_BGC_stacked_bar.png","sha256":"234226b825e3065c729197a8fe82a64a40f4732c9e65383349efa510fba3e8ab","preview_sha256":"1c0be597e35e919ee6877d069d4b86f66e539ea3fdc18f936aed70865f2db4f5","preview_content_type":"image/jpeg"}]


import matplotlib as mpl, matplotlib.pyplot as plt
print(AS.category.value_counts().to_dict())
apply_figure_style(sizes=(8,7,6))
def cls(c):
    cs=c.split(";")
    if len(cs)>1: return "Hybrid (multiple classes)"
    return {"NRPS":"NRPS","PKS":"Polyketide (PKS)","RiPP":"RiPP","terpene":"Terpene"}.get(cs[0],"Other")
AS["class"]=AS.category.apply(cls)
order_c=["NRPS","Polyketide (PKS)","RiPP","Terpene","Hybrid (multiple classes)","Other"]
cols_c=["#4C72B0","#DD8452","#55A868","#C9A227","#8172B3","#9E9E9E"]
iso_lab={"BF7":"BF7\nRhodococcoides fascians","BF3":"BF3\nKineococcus sp.","BF6":"BF6\nPsychromicrobium sp.","BF5":"BF5\nArthrobacter sp.","BF4":"BF4\nCryobacterium sp.","BF2":"BF2\nXanthobacteraceae","BF1":"BF1\nSphingomonas paucimobilis"}
ct=AS.pivot_table(index="isolate",columns="class",values="region",aggfunc="count",fill_value=0).reindex(order_iso)[[c for c in order_c if c in AS["class"].unique()]]
AS["sim"]=AS.KCB_similarity_class
order_s=["high (>=75%)","medium (50-74%)","low (<50%)","no MIBiG hit"]; cols_s=["#1B4F72","#5DADE2","#D6EAF8","#FFFFFF"]
cs=AS.pivot_table(index="isolate",columns="sim",values="region",aggfunc="count",fill_value=0).reindex(order_iso).reindex(columns=order_s,fill_value=0)
fig,axs=plt.subplots(1,2,figsize=(7.2,3.4),sharey=True,gridspec_kw={"wspace":0.08})
y=np.arange(len(order_iso))[::-1]
for ax,tab,cols_,title in [(axs[0],ct,cols_c,"BGC class"),(axs[1],cs,cols_s,"Similarity to nearest known (MIBiG) cluster")]:
    left=np.zeros(len(tab))
    for c,col in zip(tab.columns,cols_):
        ax.barh(y,tab[c].values,left=left,color=col,edgecolor="#333333",linewidth=0.4,height=0.7,label=c.replace(">=","≥"))
        left+=tab[c].values
    ax.set_title(title,loc="left"); ax.set_xlabel("Number of BGC regions")
    ax.set_xlim(0,24.5); ax.margins(y=0.02)
    ax.legend(frameon=False,loc="lower right",handlelength=1,borderaxespad=0.2)
axs[0].set_yticks(y); axs[0].set_yticklabels([iso_lab[i] for i in order_iso])
for t in axs[0].get_yticklabels():
    pass
fig.suptitle("BF7 carries the most biosynthetic gene clusters, and most regions in every isolate lack a close known match",x=0.01,ha="left",fontsize=8,y=1.02)
fig.savefig("antiSMASH_BGC_stacked_bar.png",dpi=300,bbox_inches="tight")
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
ov=[(a.get_text(),b.get_text()) for i,(a,ba) in enumerate(texts) for b,bb in texts[i+1:] if ba.overlaps(bb)]
print("overlaps:",ov[:5]); print(ct.to_string()); print(cs.to_string())
