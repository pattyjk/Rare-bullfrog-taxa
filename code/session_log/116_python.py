# session cell 116 | 2026-10-17 21:35 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_source_host_abundance.png","sha256":"4f95cb656dcaf65f8b1a512fd38ac9dfaa63c2fd19dae52f0bc9e0ecd0fef83c","preview_sha256":"fd3c7cce96bb5215432d87204d2c844693eb353c73a15f9377226cedd76792e3","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_004_AmphiBac_source_host_abundance.png","sha256":"5f539b9a644b3917c1224feec9f72724a8c6f4656bf9c80707cc4208bfee2d01","preview_sha256":"5f539b9a644b3917c1224feec9f72724a8c6f4656bf9c80707cc4208bfee2d01","preview_content_type":"image/png"}]


a2=axs[1]
for ln in a2.lines: ln.set_ydata([0,0.72])
for t in list(a2.texts):
    if t.get_text().startswith("taxon"): t.remove()
a2.text(12,215,f"taxon not\ndetected on\nsource host\n(n={int(z0.sum())})",fontsize=6,va="center")
a2.set_ylim(0,440)
a2.legend(frameon=False,loc="upper right",fontsize=6,bbox_to_anchor=(1.02,1.0))
fig.savefig("AmphiBac_source_host_abundance.png",dpi=300,bbox_inches="tight")
host.view_image("AmphiBac_source_host_abundance.png",crop=(0.5,0,1,1))
