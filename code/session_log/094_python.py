# session cell 94 | 2026-10-17 17:53 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_rank_abundance.png","sha256":"4ef4c0ae2539ab4eb2b6c81a8a674600574e05aaf4a78be7999c7bb48489a1fe","preview_sha256":"8e990a1b564be784ae8e5b8698a2321306b816a765a202120a7ad8a2ff3cfd56","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_002_AmphiBac_rank_abundance.png","sha256":"995cadd424cfdb882528b52f8502597f1b63e41a19a2c6769ffe454df251e20e","preview_sha256":"995cadd424cfdb882528b52f8502597f1b63e41a19a2c6769ffe454df251e20e","preview_content_type":"image/png"}]


for coll in ax.collections:
    if coll.get_label().startswith("AmphiBac inhibitory"): coll.set_zorder(3)
ax.set_title("Cultured isolates mostly match top-ranked skin sOTUs",loc="left")
ax.text(0.98,0.97,"BF7, BF2: V4 not found\nin meta-analysis data",transform=ax.transAxes,ha="right",va="top",fontsize=6)
fig.savefig("AmphiBac_rank_abundance.png",dpi=300,bbox_inches="tight")
host.view_image("AmphiBac_rank_abundance.png",crop=(0,0,0.66,1))
