# session cell 79 | 2026-10-17 17:48 | env: python | language: python
# files written: []


m=pd.read_csv("amb/meta.tsv",sep="\t",dtype=str,low_memory=False)
print(m["Bd.Inhibition"].value_counts(dropna=False).to_dict())
print(pd.crosstab(m["AntiFungal.Function.with.AP.function"],m["Tested.Against.Bd"].fillna("NA")))
print(m["Fungi.Tested"].value_counts().head(6).to_dict())
fa={};h=None
for l in open("amb/seqs.fasta"):
    l=l.strip()
    if l.startswith(">"): h=l[1:]; fa[h]=[]
    elif h: fa[h].append(l.upper())
fa={k:"".join(v).replace("-","") for k,v in fa.items()}
L=pd.Series({k:len(v) for k,v in fa.items()}); print(L.describe().round(0).to_dict())
print("FastaID in fasta:",m.FastaID.isin(list(fa)).sum())
# index sOTUs
sotus=list(tax.index[tax.index.str.len()==90])
K=24; idx=collections.defaultdict(list)
for s in sotus: idx[s[:K]].append(s)
def find(seq):
    hits=set()
    for qq in (seq, seq.translate(rc)[::-1]):
        for i in range(len(qq)-89):
            lst=idx.get(qq[i:i+K])
            if lst:
                w=qq[i:i+90]
                for s in lst:
                    if s==w: hits.add(s)
        if hits: break
    return hits
import time; t0=time.time()
match={k:find(v) for k,v in fa.items()}
nm=pd.Series({k:len(v) for k,v in match.items()})
print("time",round(time.time()-t0),"s; isolates with >=1 exact sOTU:",(nm>0).sum(),"of",len(nm), "multi:",(nm>1).sum())
# V4 coverage: does isolate contain 515F primer region?
hasF=pd.Series({k:bool(pr.search(v) or pr.search(v.translate(rc)[::-1])) for k,v in fa.items()})
print("isolates spanning 515F:",hasF.sum(), "| of those matched:",(nm[hasF]>0).sum())
