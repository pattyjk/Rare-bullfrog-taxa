# session cell 80 | 2026-10-17 17:48 | env: python | language: python
# files written: []


def v4win(seq):
    for qq in (seq, seq.translate(rc)[::-1]):
        mm=pr.search(qq)
        if mm and len(qq)>=mm.end()+90: return qq[mm.end():mm.end()+90]
    return None
win={k:v4win(v) for k,v in fa.items()}
res={}
for k,hs in match.items():
    if not hs: res[k]=None
    elif len(hs)==1: res[k]=next(iter(hs))
    else: res[k]=win[k] if win[k] in hs else "AMBIG"
r=pd.Series(res)
print("multi resolved by primer window:",sum(1 for k in nm[nm>1].index if r[k] not in (None,"AMBIG")),"ambiguous:",(r=="AMBIG").sum())
# example multi positions
k=nm[nm>1].index[0]; s=fa[k]
print([ (s.find(h), s.translate(rc)[::-1].find(h)) for h in match[k]], win[k] in match[k])
norm={"Inhibitory":"Inhibitory","inhibitory":"Inhibitory","Inhibitor":"Inhibitory","NonInhibitory":"Non-inhibitory","NoEffect":"Non-inhibitory","Non-inhibitor":"Non-inhibitory","No Effect":"Non-inhibitory","No effect":"Non-inhibitory","ns":"Non-inhibitory","Facilitating":"Facilitating","enhancing":"Facilitating"}
m["Bd_function"]=m["Bd.Inhibition"].map(norm)
m["sOTU"]=m.FastaID.map(r)
print(m.Bd_function.value_counts(dropna=False).to_dict())
