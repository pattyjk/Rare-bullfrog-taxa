# session cell 96 | 2026-10-17 17:58 | env: python | language: python
# files written: []


one_=nd[nd.min_mismatch==1]
amb_=one_.FastaID.map(lambda k:bool(set(win[k])-set("ACGT"))); print("1-mm with ambiguity code (N/R/Y...) in V4:",int(amb_.sum()),"of",len(one_))
pos=[]; 
for k,s in zip(one_.FastaID,one_.best_sOTU):
    w=win[k]; pos+= [i for i in range(90) if w[i]!=s[i]]
pc=pd.Series(pos).value_counts().sort_index()
print("mismatch at pos 0-4:",int(pc[pc.index<5].sum()),"85-89:",int(pc[pc.index>=85].sum()),"total",int(pc.sum()))
print("per-study (Reference) top:",one_.Reference.value_counts().head(5).to_dict())
