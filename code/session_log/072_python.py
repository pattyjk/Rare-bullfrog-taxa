# session cell 72 | 2026-10-15 22:06 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/antiSMASH_BGC_stacked_bar.png","sha256":"85a08a0fe64ee9cb10e0c2d50252ff629b58e5fb30fc936b6ce716d84725ba36","preview_sha256":"bca7812b4b2a4b6308d15d5bb52ac93b9db13ed278d2b0378c864801277a5e51","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/_view_001_antiSMASH_BGC_stacked_bar.png","sha256":"dcb5ad448b76ec9ac002c3a66f9f40098117c6b23628907c9156cd10cab8e337","preview_sha256":"dcb5ad448b76ec9ac002c3a66f9f40098117c6b23628907c9156cd10cab8e337","preview_content_type":"image/png"}]


def it(s): return r"$\it{"+s.replace(" ",r"\ ")+"}$"
iso_lab={"BF7":"BF7\n"+it("Rhodococcoides fascians"),"BF3":"BF3\n"+it("Kineococcus")+" sp.","BF6":"BF6\n"+it("Psychromicrobium")+" sp.",
 "BF5":"BF5\n"+it("Arthrobacter")+" sp.","BF4":"BF4\n"+it("Cryobacterium")+" sp.","BF2":"BF2\n"+it("Xanthobacteraceae"),"BF1":"BF1\n"+it("Sphingomonas paucimobilis")}
axs[0].set_yticklabels([iso_lab[i] for i in order_iso])
fig.suptitle("BF7 carries the most BGC regions; only 3 of 57 regions closely match a known (MIBiG) cluster",x=0.01,ha="left",fontsize=8,y=1.02)
fig.savefig("antiSMASH_BGC_stacked_bar.png",dpi=300,bbox_inches="tight")
host.view_image("antiSMASH_BGC_stacked_bar.png",crop=None)
