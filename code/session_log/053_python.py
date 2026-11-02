# session cell 53 | 2026-10-15 21:23 | env: python | language: python
# files written: []


kall=set(ti.index[ti.genus_tax=="Kineococcus"]); k97=set(ti.index[(ti.genus_tax=="Kineococcus")&(ti.Per_ident>=97)])
kex=set(s for s in set(k90.values()) if s in seqset); kexu=set(s for s in kex if v2g[s]=={"Kineococcus"})
need=(kall|kex)-set(cnt.index)
extra={}
with open(fn) as f:
    f.readline()
    for line in f:
        i=line.find("\t")
        if line[:i] in need:
            fl=line.rstrip("\r\n").split("\t")
            if len(fl)==ns+1: extra[fl[0]]=np.array(fl[1:],dtype=np.float32)
if extra: cnt=pd.concat([cnt,pd.DataFrame.from_dict(extra,orient="index",columns=samples)])
present=cnt.index
rel=cnt[good].div(T[good],axis=1)
print("Kineococcus sOTUs: tax any",len(kall),"tax>=97",len(k97),"exactV4",len(kex),"| in amphib table:",len((kall|kex)&set(present)))
Kmeth={"BLAST top hit, >=97% id":k97,"BLAST top hit, any id":kall,"Exact V4 match to NCBI 16S (any)":kex,"Exact V4 match, unique to genus":kexu}
for lab,ids in Kmeth.items():
    R=pd.concat([R[~((R.genus=="Kineococcus")&(R.method==lab))],pd.DataFrame([dict(genus="Kineococcus",method=lab,**summ_f) for summ_f in [summ.__class__ and None] if False])])
newrows=[dict(genus="Kineococcus",method=lab,**(lambda ids:(summ_amph(ids)))( [s for s in ids if s in present])) for lab,ids in []]
