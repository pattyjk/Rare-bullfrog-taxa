# session cell 105 | 2026-10-17 18:13 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_cultured_vs_abundance.png","sha256":"3b2fa904626b71ff8782cd54017548e9f63d5bba824871a751965b1733ca7b20","preview_sha256":"1b2441319296bf46396f77288469975ada6bdfc6e0c397b24038019281801df3","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_008_AmphiBac_cultured_vs_abundance.png","sha256":"53478ea79271409f5e655ce3d7ffc6f25386eda17a04f6f5986519c04ef61256","preview_sha256":"53478ea79271409f5e655ce3d7ffc6f25386eda17a04f6f5986519c04ef61256","preview_content_type":"image/png"}]


ax3.get_legend().remove()
for t in list(ax3.texts): t.remove()
ax3.text(0.97,0.04,f"Line: median per 10-sOTU bin\nSpearman ρ = {rho:.2f}, n = {len(xm):,}",transform=ax3.transAxes,ha="right",va="bottom",fontsize=6,bbox=dict(facecolor="white",edgecolor="none",alpha=0.8,pad=1))
ax3.set_title("Read share declines above ~20 matches",loc="left")
fig.savefig("AmphiBac_cultured_vs_abundance.png",dpi=300,bbox_inches="tight")
host.view_image("AmphiBac_cultured_vs_abundance.png",crop=(0.64,0,1,1))
print([ (lo, int(((xm>=lo)&(xm<lo+10)).sum())) for lo in bins[:3]])
