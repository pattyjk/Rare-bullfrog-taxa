# session cell 123 | 2026-10-18 11:30 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_003_AmphiBac_abundance_class_function.png","sha256":"d2d1296e8e4d07095067b840b410bf0969150b6dd0a1fedc9a499235282bdf2b","preview_sha256":"d2d1296e8e4d07095067b840b410bf0969150b6dd0a1fedc9a499235282bdf2b","preview_content_type":"image/png"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_abundance_class_function.png","sha256":"7e273c662858e5d5c54aeb328e96b99f8ac5268d66e74d2d5bf5770440bd7e89","preview_sha256":"0fe5493ed64e67f4434d5538a7d975e514bbe913a24084acfb22f8c1a6098e7c","preview_content_type":"image/jpeg"}]


for t in list(axs[1].texts):
    if t.get_text().startswith("χ²"): t.remove()
axs[1].text(0.5,0.855,f"inhibitory vs not: χ² p = {p_chi:.3f}",transform=axs[1].transAxes,ha="center",fontsize=6)
fig.savefig("AmphiBac_abundance_class_function.png",dpi=300,bbox_inches="tight")
host.view_image("AmphiBac_abundance_class_function.png",crop=(0.33,0,0.67,0.4))
