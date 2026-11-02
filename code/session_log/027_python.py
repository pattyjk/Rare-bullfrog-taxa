# session cell 27 | 2026-10-15 19:59 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/sub_counts.pkl","sha256":"f1bf79cb37c7a2bc13780216e546279fabf370d8cecc17f12f2f6390a667ddcc"}]


need=set(exm.sOTU)-set(sub.index)
extra={}
with open(fn) as f:
    f.readline()
    for line in f:
        i=line.find("\t")
        if line[:i] in need:
            fl=line.rstrip("\r\n").split("\t")
            if len(fl)==ns+1: extra[fl[0]]=np.array(fl[1:],dtype=np.float32)
print(len(need),len(extra))
cnt=pd.concat([sub,pd.DataFrame.from_dict(extra,orient="index",columns=samples)])
cnt.to_pickle("meta/sub_counts.pkl")
