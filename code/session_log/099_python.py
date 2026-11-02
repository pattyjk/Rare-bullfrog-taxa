# session cell 99 | 2026-10-17 17:59 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_cultured_vs_abundance.png","sha256":"52d3e4a12829ed7158cf71296be9099f4c8f4bdc7869c4c08491431fb40a86ec","preview_sha256":"a324cd7da7c8d82b3c6fee777f0439ef2931f323b39256d6b45494cecf7b7568","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_004_AmphiBac_cultured_vs_abundance.png","sha256":"785718deccba85d3fcafc75a60767c6150b9117361d5b61bed31ebd04082306c","preview_sha256":"785718deccba85d3fcafc75a60767c6150b9117361d5b61bed31ebd04082306c","preview_content_type":"image/png"}]


axs[0].set_title("All dominant sOTUs have been cultured",loc="left")
axs[1].set_title("Few sOTUs, most reads",loc="left")
for t in list(axs[0].texts):
    if t.get_text() in ("428/48,334","476"): t.remove()
axs[0].text(2-w/2-0.02,4,"428/48.3k\n(0.9%)",ha="right",va="bottom",fontsize=6)
axs[0].text(2+w/2+0.02,4,"476\n(1.0%)",ha="left",va="bottom",fontsize=6)
fig.savefig("AmphiBac_cultured_vs_abundance.png",dpi=300,bbox_inches="tight")
r=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print([(a_.get_text(),b_.get_text()) for i,(a_,ba) in enumerate(texts) for b_,bb in texts[i+1:] if ba.overlaps(bb)][:8])
host.view_image("AmphiBac_cultured_vs_abundance.png")
