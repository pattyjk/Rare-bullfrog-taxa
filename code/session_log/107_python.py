# session cell 107 | 2026-10-17 21:32 | env: python | language: python
# files written: []


import pandas as pd, numpy as np, re, collections
from scipy import stats as sst
rc=str.maketrans("ACGT","TGCA"); pr=re.compile("GTG[CT]CAGC[AC]GCCGCGGTAA")
ST=pd.read_pickle("meta/sOTU_stats.pkl")
presS=list(ST.index[ST.n_samples>0])
m=pd.read_csv("amb/meta.tsv",sep="\t",dtype=str,low_memory=False)
norm={"Inhibitory":"Inhibitory","inhibitory":"Inhibitory","Inhibitor":"Inhibitory","NonInhibitory":"Non-inhibitory","NoEffect":"Non-inhibitory","Non-inhibitor":"Non-inhibitory","No Effect":"Non-inhibitory","No effect":"Non-inhibitory","ns":"Non-inhibitory","Facilitating":"Facilitating","enhancing":"Facilitating"}
m["Bd_function"]=m["Bd.Inhibition"].map(norm)
fa={};h=None
for l in open("amb/seqs.fasta"):
    l=l.strip()
    if l.startswith(">"): h=l[1:]; fa[h]=[]
    elif h: fa[h].append(l.upper())
fa={k:"".join(v).replace("-","") for k,v in fa.items()}
K=24; idx2=collections.defaultdict(list)
for s in presS: idx2[s[:K]].append(s)
def find2(seq):
    hits=set()
    for qq in (seq, seq.translate(rc)[::-1]):
        for i in range(len(qq)-89):
            lst=idx2.get(qq[i:i+K])
            if lst:
                w=qq[i:i+90]; hits.update(s for s in lst if s==w)
    return hits
def v4win(seq):
    for qq in (seq, seq.translate(rc)[::-1]):
        mm=pr.search(qq)
        if mm and len(qq)>=mm.end()+90: return qq[mm.end():mm.end()+90]
    return None
match={k:find2(v) for k,v in fa.items()}; win={k:v4win(v) for k,v in fa.items()}
PS=np.frombuffer("".join(presS).encode(),dtype=np.uint8).reshape(len(presS),90)
meanv=ST.loc[presS,"mean_rel_pct"].values
for k in fa:
    if not match[k] and win[k] and set(win[k])<=set("ACGTN"):
        d=(PS!=np.frombuffer(win[k].encode(),dtype=np.uint8)).sum(1); dm=d.min()
        if dm<=3:
            ii=np.where(d==dm)[0]; match[k]={presS[ii[np.argmax(meanv[ii])]]}
m["sOTUs"]=m.FastaID.map(lambda k:";".join(sorted(match[k])))
print("isolates matched (<=3mm):",(m.sOTUs!="").sum())
mp=pd.read_csv(host.artifact_path("e39afe94-4d38-4b9a-bb95-7731a3f8da7e"),sep="\t",dtype=str).set_index("#SampleID")
RM0=pd.read_pickle("meta/matched_rel.pkl"); gsamples=list(RM0.columns)
mpg=mp.loc[gsamples]
mpg["hsp"]=(mpg.host_genus.str.lower().str.strip()+" "+mpg.host_species.str.lower().str.strip())
print("meta host species:",mpg.hsp.nunique()); 
def nh(s): 
    s=str(s).replace("."," ").replace("_"," ").strip().lower(); s=re.sub(r"\s+"," ",s); return " ".join(s.split()[:2])
m["hsp"]=m["Host.Species"].map(nh)
ab_h=m[m.sOTUs!=""].hsp.value_counts()
ov=ab_h[ab_h.index.isin(mpg.hsp.unique())]
print("AmphiBac matched host species:",len(ab_h),"| exact species in meta:",len(ov),"isolates",int(ov.sum()),"of",int(ab_h.sum()))
print("not found (top):",ab_h[~ab_h.index.isin(mpg.hsp.unique())].head(15).to_dict())
print("meta hosts sample counts top:",mpg.hsp.value_counts().head(12).to_dict())
