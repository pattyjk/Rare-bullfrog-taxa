# session cell 10 | 2026-10-15 19:48 | env: python | language: python
# files written: []


import pandas as pd, re, numpy as np
tax=pd.read_csv("meta/tax.txt",sep="\t",dtype=str)
print(tax.shape, tax.sOTU.nunique(), tax.sOTU.str.len().describe()[["min","50%","max"]].tolist())
tax["Per_ident"]=tax.Per_ident.astype(float)
tax["genus_tax"]=tax.Taxonomy.fillna("").str.rstrip("; ").str.split("; ").str[-1]
tax["hit_genus"]=tax.Hit_name.fillna("").str.split().str[0]
tax["hit_epithet"]=tax.Hit_name.fillna("").str.split().str[1]
rhodoides={"fascians","corynebacterioides","kroppenstedtii","kyotonensis","kyotonense","navarretei","trifolii","yunnanensis","yunnanense","luteus"}
targets=["Rhodococcoides","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter","Sphingomonas"]
def assign(r):
    lin=r.Taxonomy if isinstance(r.Taxonomy,str) else ""
    if r.hit_genus in ("Rhodococcus","Rhodococcoides") and r.hit_epithet in rhodoides: return "Rhodococcoides"
    if r.hit_genus=="Rhodococcoides" or "; Rhodococcoides" in lin: return "Rhodococcoides"
    for g in targets[1:]:
        if f"; {g};" in lin+";" or r.hit_genus==g: return g
    if "; Rhodococcus;" in lin+";": return "Rhodococcus (other)"
    return None
tax["target"]=tax.apply(assign,axis=1)
print(tax.target.value_counts())
for g in ["Rhodococcoides","Rhodococcus (other)","Psychromicrobium","Cryobacterium"]:
    s=tax[tax.target==g]
    print(g, s.Hit_name.str.split().str[:2].str.join(" ").value_counts().head(8).to_dict(), "ident", s.Per_ident.describe()[["min","50%"]].tolist())
print(tax[tax.Hit_name.str.contains("Psychromicrobium|lacuslunae",na=False)].shape)
print(tax.Per_ident.describe())
