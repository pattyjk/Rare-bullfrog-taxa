# session cell 177 | 2026-11-02 11:12 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/work/amphib_only_otu_table.txt","sha256":null,"oversize":true}]

import pandas as pd, numpy as np, re, zipfile, io, os, time
t0=time.time()
os.makedirs('work',exist_ok=True)
z=zipfile.ZipFile(host.artifact_path('9ab57aec-cb7f-47eb-a688-a9f6e2156307')); nm=[n for n in z.namelist() if n.endswith('.txt')][0]
z.extract(nm,'work'); txt='work/'+nm
with open(txt) as f: hdr=f.readline().rstrip('\n').split('\t')
nf=len(hdr)+1
ids=[];rows=[]
with open(txt) as f:
    f.readline()
    for line in f:
        p=line.rstrip('\n').split('\t')
        if len(p)!=nf or len(p[0])!=90: continue
        ids.append(p[0]); rows.append(np.array(p[1:],dtype=np.float32))
X=np.vstack(rows); del rows
tot=X.sum(0); keep=tot>=1000
R=X[:,keep]/tot[keep]*100
mean_rel=R.mean(1); prev=(R>0).mean(1)*100
st=pd.DataFrame({'sOTU':ids,'mean_rel_pct':mean_rel,'prevalence_pct':prev})
st=st[st.mean_rel_pct>0].reset_index(drop=True)
print(keep.sum(), len(ids), len(st), round(time.time()-t0))