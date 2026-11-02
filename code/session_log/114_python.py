# session cell 114 | 2026-10-17 21:34 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/AmphiBac_source_host_abundance.png","sha256":"176e938072a27d01797b30ed0140f0e84bf0a25b7e3a01e3d7250e329d0c1be7","preview_sha256":"e90919258853ea3344b7456bec7689c9e2c49561f515adfdaa7555699d9df303","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_002_AmphiBac_source_host_abundance.png","sha256":"7d0027bc059b51e0ac5a1e1ee40054570a916a53e104b50ea27d913b7001ac43","preview_sha256":"7d0027bc059b51e0ac5a1e1ee40054570a916a53e104b50ea27d913b7001ac43","preview_content_type":"image/png"}]


z0=(sp.home_group_mean==0); print("species-level units with taxon absent on source host:",int(z0.sum()),"of",len(sp), sp[z0].function.value_counts().to_dict())
hi=sp.pct_other_groups_below_home>=90; print(">=90%:",int(hi.sum()), sp[hi].function.value_counts().to_dict())
nz=sp[~z0]; print("excluding absent: median % below:",nz.groupby("function").pct_other_groups_below_home.median().round(1).to_dict(), "MWU p",sst.mannwhitneyu(nz[nz.function=="Inhibitory"].pct_other_groups_below_home,nz[nz.function!="Inhibitory"].pct_other_groups_below_home).pvalue)
print("abs share by function:",sp.groupby("function").apply(lambda d:round(100*(d.home_group_mean==0).mean(),1)).to_dict())
# are absent-at-home units from studies not in meta? check whether home samples are from same study (doi) 
axs[0].set_title("Source host vs other hosts",loc="left",fontsize=8)
axs[1].set_title("Source host: often top-ranked or taxon absent",loc="left",fontsize=8)
axs[1].annotate(f"not detected on\nsource host (n={int(z0.sum())})",xy=(5,300),xytext=(18,330),fontsize=6,arrowprops=dict(arrowstyle="-",lw=0.5))
fig.savefig("AmphiBac_source_host_abundance.png",dpi=300,bbox_inches="tight")
host.view_image("AmphiBac_source_host_abundance.png")
