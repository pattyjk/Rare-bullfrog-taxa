# session cell 115 | 2026-10-17 21:34 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_source_host_abundance.png","sha256":"a86cc864f1c713edd4cd85943459e2d0241a0227175d2dede2691bb9064d8096","preview_sha256":"12b40d628dfc32cfef6a225cc62795013a869b80f712c9e7d6ccd1b7125f4049","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_003_AmphiBac_source_host_abundance.png","sha256":"6d09bd6f109d56d3534e8cbbb0603fc57dbdbaec8ed01ca9ac678033e3d5456b","preview_sha256":"6d09bd6f109d56d3534e8cbbb0603fc57dbdbaec8ed01ca9ac678033e3d5456b","preview_content_type":"image/png"}]


for t in list(axs[1].texts): t.remove()
for t in list(axs[0].texts):
    if t.get_text().startswith("above line"): t.remove()
axs[0].text(0.97,0.06,"below line: higher on source host",transform=axs[0].transAxes,ha="right",fontsize=6,color="#444444")
axs[1].annotate(f"taxon not detected on\nsource host (n={int(z0.sum())})",xy=(10,250),xytext=(16,235),fontsize=6,va="center")
axs[1].legend(frameon=False,loc="upper right",fontsize=6,bbox_to_anchor=(1.0,1.0))
panel_letter(axs[1],"b") if not any(t.get_text()=="b" for t in fig.texts+axs[1].texts) else None
fig.savefig("AmphiBac_source_host_abundance.png",dpi=300,bbox_inches="tight")
host.view_image("AmphiBac_source_host_abundance.png",crop=(0.5,0,1,1))
