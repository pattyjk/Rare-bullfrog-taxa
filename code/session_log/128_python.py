# session cell 128 | 2026-10-18 11:33 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_007_AmphiBac_abundance_class_function_isolates.png","sha256":"f63a3e8bdd4216309fca3ff8cc2d502b94cdb4fbde5bb1c893bca8bc152f2dcd","preview_sha256":"f63a3e8bdd4216309fca3ff8cc2d502b94cdb4fbde5bb1c893bca8bc152f2dcd","preview_content_type":"image/png"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_abundance_class_function_isolates.png","sha256":"e32bb270048c0d16aed50578aea851581bce037e02dcf90ab0f77414b20847ea","preview_sha256":"c9f0103726075fba7c08fbd2b32b22093d8d813d0dfb12b73d3db22632a8d400","preview_content_type":"image/jpeg"}]


axs[0].set_title("53% of isolates are rare taxa",loc="left")
for c_ in list(axs[2].collections): c_.remove()
lo_,hi_=-215,215
rng=np.random.default_rng(3)
for i,d_ in enumerate(di,1):
    dd=np.clip(d_,lo_+5,hi_-5); axs[2].scatter(i+rng.uniform(-0.18,0.18,len(dd)),dd,s=1.2,color="#2166AC",alpha=0.18,lw=0,zorder=3)
axs[2].set_ylim(lo_,hi_)
for t in list(axs[2].texts):
    if t.get_text().startswith("Kruskal"): t.set_position((0.98,0.99))
fig.savefig("AmphiBac_abundance_class_function_isolates.png",dpi=300,bbox_inches="tight")
host.view_image("AmphiBac_abundance_class_function_isolates.png")
print("points clipped to axis edge:",int(((allsc<lo_+5)|(allsc>hi_-5)).sum()))
