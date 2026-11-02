# session cell 125 | 2026-10-18 11:30 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_005_AmphiBac_abundance_class_function.png","sha256":"f3174a381aa1d4e7cf05259b7670d1839473259099079ab9bc286a7c1b7629a4","preview_sha256":"f3174a381aa1d4e7cf05259b7670d1839473259099079ab9bc286a7c1b7629a4","preview_content_type":"image/png"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_abundance_class_function.png","sha256":"76dd5358b6e79c24c3ff45a650b5df71c7d6a83aed48bcc5b5e3e06de165524c","preview_sha256":"988e4283f8d86dd6d3e2aa35684143b45b60729b30e8ffcb2fc663d3b6083fb2","preview_content_type":"image/jpeg"}]


h_,l_=axs[1].get_legend_handles_labels()
axs[1].legend(h_,["Inhibitory","No effect/mixed","Enhancing"],frameon=False,loc="upper center",ncol=3,fontsize=6,handlelength=0.8,handletextpad=0.3,columnspacing=0.5,borderaxespad=0,bbox_to_anchor=(0.52,1.0))
fig.savefig("AmphiBac_abundance_class_function.png",dpi=300,bbox_inches="tight")
host.view_image("AmphiBac_abundance_class_function.png",crop=(0.33,0,0.7,0.25))
