# session cell 81 | 2026-10-17 17:49 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/sOTU_stats.pkl","sha256":"ea4c266627bdc8a19049f47b6f534851ac194618d2ab9ea82af1320094326e18"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/matched_rel.pkl","sha256":"d8fc030b1b0c263462d1435d29ded69190b53a78173dea7ece47cd32e2b117c2"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/otu/amphib_only_otu_table.txt","sha256":null,"oversize":true}]


import zipfile, io, time; t0=time.time()
zp=host.artifact_path("9ab57aec-cb7f-47eb-a688-a9f6e2156307")
zipfile.ZipFile(zp).extractall("meta/otu")
fn="meta/otu/amphib_only_otu_table.txt"
with open(fn) as f: samples=f.readline().rstrip("\n").split("\t")
ns=len(samples); tot=np.zeros(ns)
with open(fn) as f:
    f.readline()
    for line in f:
        fl=line.rstrip("\r\n").split("\t")
        if len(fl)!=ns+1 or len(fl[0])!=90: continue
        try: tot+=np.array(fl[1:],dtype=np.float64)
        except ValueError: pass
good=tot>=1000; ng=good.sum()
need=set().union(*[h for h in match.values() if h])
stats={}; keepm={}
with open(fn) as f:
    f.readline()
    for line in f:
        fl=line.rstrip("\r\n").split("\t")
        if len(fl)!=ns+1 or len(fl[0])!=90: continue
        try: v=np.array(fl[1:],dtype=np.float64)
        except ValueError: continue
        rr=v[good]/tot[good]
        stats[fl[0]]=(v[good].sum(), (rr>0).sum(), rr.mean(), rr.max())
        if fl[0] in need: keepm[fl[0]]=rr.astype(np.float32)
ST=pd.DataFrame.from_dict(stats,orient="index",columns=["reads","n_samples","mean_rel","max_rel"])
ST["prevalence_pct"]=100*ST.n_samples/ng; ST["mean_rel_pct"]=100*ST.mean_rel; ST["max_rel_pct"]=100*ST.max_rel
gsamples=[s for s,g in zip(samples,good) if g]
RM=pd.DataFrame.from_dict(keepm,orient="index",columns=gsamples)*100
ST.to_pickle("meta/sOTU_stats.pkl"); RM.to_pickle("meta/matched_rel.pkl")
print("samples used",ng,"sOTUs",len(ST),"matched sOTUs present",len(RM),"of",len(need), round(time.time()-t0),"s")
