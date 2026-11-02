# session cell 110 | 2026-10-17 21:33 | env: python | language: python
# files written: []


tot=np.zeros(ns); raw={}
with open(fn) as f:
    f.readline()
    for line in f:
        fl=line.rstrip("\r\n").split("\t")
        if len(fl)!=ns+1 or len(fl[0])!=90: continue
        try: v=np.array(fl[1:],dtype=np.float64)
        except ValueError: continue
        tot+=v
        if fl[0] in miss: raw[fl[0]]=v
ci=np.array(colidx)
add=pd.DataFrame({k:(v[ci]/tot[ci]*100).astype(np.float32) for k,v in raw.items()},index=gsamples).T
RM=pd.concat([RM0,add]); RM=RM.loc[~RM.index.duplicated()]
print(RM.shape, all(s in RM.index for s in need))
assert {'RM','mg','mpg'} <= globals().keys()
