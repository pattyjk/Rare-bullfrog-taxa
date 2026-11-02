# session cell 118 | 2026-10-17 21:35 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_source_host_abundance.png","sha256":"aa313aa1296db5f6cccef862d1ea2067f2e96c6814823aa1c2866c9e9664afc8","preview_sha256":"6f174ec6662b9a4fca15913fd3a0ba08c490b0c24369152764ccd5bdbdac73c9","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_006_AmphiBac_source_host_abundance.png","sha256":"bb2d32ec60f7994352589425660363ca51c758927714ee48547bd9a8edeba3da","preview_sha256":"bb2d32ec60f7994352589425660363ca51c758927714ee48547bd9a8edeba3da","preview_content_type":"image/png"}]


h_,l_=a2.get_legend_handles_labels()
a2.legend(h_,["Inhibitory (n=684; median 58%)","Non-inhibitory/facil. (n=792; 50%)"],frameon=False,loc="upper right",fontsize=6,bbox_to_anchor=(1.0,1.0))
fig.savefig("AmphiBac_source_host_abundance.png",dpi=300,bbox_inches="tight")
host.view_image("AmphiBac_source_host_abundance.png",crop=(0.5,0,1,0.3))
