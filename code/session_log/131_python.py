# session cell 131 | 2026-10-18 11:37 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/RMh.pkl","sha256":"95b4cae05523ea21a2e1544deb0175033149b1b492f7bc295634dd39e7c25e98"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/otu/amphib_only_otu_table.txt","sha256":null,"oversize":true}]


import zipfile, os
os.makedirs("meta/otu",exist_ok=True)
zipfile.ZipFile(host.artifact_path("9ab57aec-cb7f-47eb-a688-a9f6e2156307")).extractall("meta/otu")
fn="meta/otu/amphib_only_otu_table.txt"
with open(fn) as f: samples=f.readline().rstrip("\n").split("\t")
ns=len(samples); need=set(T.index); tot=np.zeros(ns); raw={}
with open(fn) as f:
    f.readline()
    for line in f:
        fl=line.rstrip("\r\n").split("\t")
        if len(fl)!=ns+1 or len(fl[0])!=90: continue
        try: v=np.array(fl[1:],dtype=np.float64)
        except ValueError: continue
        tot+=v
        if fl[0] in need: raw[fl[0]]=v
good=tot>=1000; gs=[s for s,g in zip(samples,good) if g]
RMh=pd.DataFrame({k:(v[good]/tot[good]*100) for k,v in raw.items()},index=gs).T
RMh.to_pickle("meta/RMh.pkl")
print(RMh.shape, len(need))
