# session cell 26 | 2026-10-15 19:59 | env: python | language: python
# files written: []


import json
refs=json.load(open("meta/ref/v4refs.json"))
v2g={}
for g,d in refs.items():
    for h,s in d.items(): v2g.setdefault(s,set()).add(g)
seqset=set(tax.sOTU)
ex={s:gs for s,gs in v2g.items() if s in seqset}
tgt6=["Rhodococcoides","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter","Sphingomonas"]
rows=[]
for s,gs in ex.items():
    t=tax[tax.sOTU==s].iloc[0]
    rows.append(dict(sOTU=s,ref_genera=";".join(sorted(gs)),blast_hit=" ".join(str(t.Hit_name).split()[:2]),blast_ident=t.Per_ident,blast_target=t.target,in_sub=s in sub.index))
exm=pd.DataFrame(rows)
for g in tgt6:
    e=exm[exm.ref_genera.str.contains(g)]
    print(g,"exact sOTUs:",len(e),"unique-to-genus:",(e.ref_genera==g).sum(),"| blast top hits:",e.blast_hit.value_counts().head(4).to_dict(),"| all in sub:",e.in_sub.all())
