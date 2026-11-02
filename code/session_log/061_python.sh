# session cell 61 | 2026-10-15 21:32 | env: python | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF2.min1k.fasta","sha256":"fe90c529576c87b5204d2e3498066dca1c94d95c4527a2e90fe69ac23caca8e4"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF1.min1k.fasta","sha256":"acca257683e54a80c832a5b1e5fc8f4ae7e308e3c0bcb154582066fe3c4db139"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF5.min1k.fasta","sha256":"86d940600df56bc75ea6c46bc7a5b38362904bc587d77ebe057da3f30a9d50be"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF3.min1k.fasta","sha256":"732790f58d0f9db374a4e1be9ec7d313848dd24fa976be82f647f1ec10b1534e"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF4.min1k.fasta","sha256":"1ddca234640aa5e435722ed8bda87ad7d47086778f7220f257e831787404cd96"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF7.min1k.fasta","sha256":"c12a84314f448db532f06ca041774a294948c360742e5bdb336e10fdfd133487"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF6.min1k.fasta","sha256":"1267e8fa72f6d3506e17408f560075be572b329efef4f9c9fc26adcc009289ea"}]

cd as; python - <<'E'
import glob
for f in sorted(glob.glob("*.fasta")):
    if f.endswith(".min1k.fasta"): continue
    recs=[];h=None
    for l in open(f):
        l=l.rstrip()
        if l.startswith(">"): h=l[1:]; recs.append([h,[]])
        else: recs[-1][1].append(l)
    recs=[(h,"".join(s)) for h,s in recs]
    k=[(h,s) for h,s in recs if len(s)>=1000]
    with open(f.replace(".fasta",".min1k.fasta"),"w") as o:
        for h,s in k: o.write(f">{h}\n"+"\n".join(s[i:i+80] for i in range(0,len(s),80))+"\n")
    print(f, len(k), sum(len(s) for _,s in k))
E