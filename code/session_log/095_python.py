# session cell 95 | 2026-10-17 17:58 | env: python | language: python
# files written: []


nd=m[m.detection.str.startswith("not detected")].copy()
PS=np.frombuffer("".join(presS).encode(),dtype=np.uint8).reshape(len(presS),90)
best_d=[];best_s=[];n_best=[]
for k in nd.FastaID:
    w=win[k]; a_=np.frombuffer(w.encode(),dtype=np.uint8)
    d=(PS!=a_).sum(1); dm=int(d.min()); ii=np.where(d==dm)[0]
    best_d.append(dm); n_best.append(len(ii)); best_s.append(presS[ii[np.argmax(ST.mean_rel_pct.values[[ST.index.get_loc(presS[j]) for j in ii]])]])
nd["min_mismatch"]=best_d; nd["best_sOTU"]=best_s; nd["n_tied_sOTUs"]=n_best
nd["pct_identity"]=(100*(90-nd.min_mismatch)/90).round(1)
print(nd.min_mismatch.value_counts().sort_index().head(12).to_dict())
cum={t:int((nd.min_mismatch<=t).sum()) for t in [1,2,3,4,5]}; print("cumulative rescued:",cum)
# non-V4 isolates: what do they cover?
unk=m[m.detection.str.startswith("unknown")]
L_unk=unk.FastaID.map(lambda k:len(fa[k])); print("unknown: length median",L_unk.median())
# Does the ≥99% (<=1 mm) rescue change the inhibition comparison?
nd["rescued_1mm"]=nd.min_mismatch<=1
b1=nd[nd.Bd_function.notna()]
print(b1.groupby("rescued_1mm").Bd_function.apply(lambda x:round((x=="Inhibitory").mean(),3)).to_dict(), b1.rescued_1mm.value_counts().to_dict())
det1=bt.detection.str.startswith("detected")|bt.FastaID.isin(nd.FastaID[nd.rescued_1mm])
kn1=bt[~bt.detection.str.startswith("unknown")]
d1=kn1.FastaID.isin(bt.FastaID[det1])
print("with <=1mm counted as detected: inhib frac detected",round((kn1[d1].Bd_function=="Inhibitory").mean(),3),"not",round((kn1[~d1].Bd_function=="Inhibitory").mean(),3), sst.chi2_contingency(pd.crosstab(d1,kn1.Bd_function=="Inhibitory"))[1])
print(nd[nd.rescued_1mm].best_sOTU.map(ST.mean_rel_pct).describe().round(4).to_dict())
