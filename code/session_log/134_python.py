# session cell 134 | 2026-10-18 11:38 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_abundance_class_definitions.png","sha256":"97027490c2063590eadbb53ecba2c64f87f111106ba299abf1145fb1198dcb3a","preview_sha256":"2bc54869420a21efadb7bf668aa044d1a95691a101c93345491b28f615117846","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_009_AmphiBac_abundance_class_definitions.png","sha256":"80bd1fd40c90886d434191e0d94b7b32362c7ddd4ce172d46d9daa3e1cc40f50","preview_sha256":"80bd1fd40c90886d434191e0d94b7b32362c7ddd4ce172d46d9daa3e1cc40f50","preview_content_type":"image/png"}]


axs[0,1].get_legend().remove()
h_,l_=axs[0,1].get_legend_handles_labels()
fig.legend(h_,l_,frameon=False,loc="lower center",ncol=3,fontsize=7,bbox_to_anchor=(0.5,0.985),handlelength=0.9)
for r in range(2): axs[r,1].set_ylabel("Taxa by Bd\nfunction (%)")
axs[2,1].set_ylabel("Units by Bd\nfunction (%)")
lab4=["Dom.\n≥1%","Interm.\n0.1–1%","Rare\n<0.1%","Not\ndet."]
for c in range(3): axs[2,c].set_xticklabels(lab4)
for r in range(3):
    axs[r,2].set_ylim(-200,185)
    for t in axs[r,2].texts:
        if t.get_text().startswith("KW"): t.set_position((0.98,0.99))
fig.text(0.5,-0.005,"Not det. = taxon not detected in any sample of its source host species. Numbers under boxes = units with an inhibition score (96-well CFS assays only).",ha="center",fontsize=6,color="#444444")
fig.savefig("AmphiBac_abundance_class_definitions.png",dpi=300,bbox_inches="tight")
r_=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r_)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print([(a_.get_text()[:25],b_.get_text()[:25]) for i,(a_,ba) in enumerate(texts) for b_,bb in texts[i+1:] if ba.overlaps(bb) and len(a_.get_text())>1 and len(b_.get_text())>1][:12])
host.view_image("AmphiBac_abundance_class_definitions.png")
