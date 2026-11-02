# session cell 75 | 2026-10-17 16:58 | env: python | language: python
# files written: []


import pandas as pd, numpy as np, re, collections, json
rc=str.maketrans("ACGT","TGCA"); pr=re.compile("GTG[CT]CAGC[AC]GCCGCGGTAA")
tgt=["Rhodococcoides","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter","Sphingomonas","Kineococcus","Rhodococcus (other)"]
rhodoides={"fascians","corynebacterioides","kroppenstedtii","kyotonensis","kyotonense","navarretei","trifolii","yunnanensis","yunnanense","luteus"}
gens=["Rhodococcoides","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter","Sphingomonas","Kineococcus","Rhodococcus","Paenarthrobacter","Pseudarthrobacter","Glutamicibacter"]
# reference 120-bp variants
ref120={}
for g in gens:
    seqs={};h=None
    for l in open(f"meta/ref/{g}.fasta"):
        l=l.strip()
        if l.startswith(">"): h=l[1:]; seqs[h]=[]
        elif h: seqs[h].append(l.upper())
    for h,s in seqs.items():
        s="".join(s)
        for qq in (s,s.translate(rc)[::-1]):
            mm=pr.search(qq)
            if mm:
                if len(qq)>=mm.end()+120:
                    fr=qq[mm.end():mm.end()+120]
                    if set(fr)<=set("ACGT"): ref120.setdefault(fr,{}).setdefault(g,collections.Counter())[" ".join(h.split()[1:3])]+=1
                break
lab=lambda g:"Rhodococcus (other)" if g=="Rhodococcus" else g
# ASVs
q={};h=None
for l in open("sal/asvs.fasta"):
    l=l.strip()
    if l.startswith(">"): h=l[1:]; q[h]=""
    else: q[h]+=l.upper()
tab=pd.read_csv("sal/asv_table.txt",sep="\t",skiprows=1,index_col=0)
print(tab.shape, tab.index.isin(list(q)).sum())
Tt=tab.sum(); print("reads/sample:",Tt.describe().round(0).to_dict())
RB=tab.div(Tt,axis=1)*100
# BLAST
b=pd.read_csv("sal/asv_blast.tsv",sep="\t",header=None,names="qseqid sacc pident length mismatch gapopen qstart qend evalue bitscore staxids sscinames stitle".split())
def gen_of(n):
    w=str(n).split()
    if len(w)>=2 and w[0] in("Rhodococcus","Rhodococcoides") and w[1] in rhodoides: return "Rhodococcoides"
    if w and w[0]=="Rhodococcus": return "Rhodococcus (other)"
    return w[0] if w else None
b["genus"]=b.sscinames.apply(gen_of)
best=b[b.bitscore==b.groupby("qseqid").bitscore.transform("max")]
first=b.groupby("qseqid",sort=False).head(1).set_index("qseqid")
B=first[["sscinames","genus","pident","length"]].rename(columns={"sscinames":"blast_top_hit","genus":"blast_top_genus","pident":"blast_pct_identity","length":"aln_length"})
B=B.join(best.groupby("qseqid").agg(tied_genera=("genus",lambda x:";".join(sorted(set(map(str,x))))),n_tied_hits=("sacc","size"),tied_species=("sscinames",lambda x:"; ".join(pd.Series(x).value_counts().index[:5]))))
C=B.reindex(tab.index)
C["sequence"]=[q[a] for a in tab.index]
C["exact120_ref_genera"]=C.sequence.apply(lambda s:";".join(sorted(lab(g) for g in ref120.get(s,{}))))
C["exact120_n_refs"]=C.sequence.apply(lambda s:sum(sum(c.values()) for c in ref120.get(s,{}).values()))
C["exact120_ref_species"]=C.sequence.apply(lambda s:"; ".join(f"{sp} ({n})" for g,c in ref120.get(s,{}).items() for sp,n in c.most_common(3)))
# meta-analysis taxonomy (90-bp prefix)
tax=pd.read_csv(host.artifact_path("9439c2bb-3d6a-4202-bb92-68b1ad9c937c"),sep="\t",dtype=str).drop_duplicates("sOTU").set_index("sOTU")
tax["Per_ident"]=tax.Per_ident.astype(float)
tax["genus_tax"]=tax.Taxonomy.fillna("").str.rstrip("; ").str.split("; ").str[-1]
def meta_gen(s):
    if s not in tax.index: return ""
    t=tax.loc[s]; w=str(t.Hit_name).split()
    if t.Per_ident<97: return ""
    if len(w)>1 and w[0]=="Rhodococcus" and w[1] in rhodoides: return "Rhodococcoides"
    g=t.genus_tax
    return "Rhodococcus (other)" if g=="Rhodococcus" else g
p90=C.sequence.str[:90]
C["metaTaxonomy90"]=p90.apply(meta_gen)
C["metaTaxonomy90_hit"]=p90.apply(lambda s:f"{str(tax.loc[s].Hit_name).split(' 16S')[0]} ({tax.loc[s].Per_ident}%)" if s in tax.index else "")
# amphibian prevalence from saved sOTU_detail sheet
amv=host.artifacts(filename="amphibian_genus_abundance.xlsx",exact=True)["artifacts"][0]["latest_version_id"]
sd=pd.read_excel(host.artifact_path(amv),sheet_name="sOTU_detail").drop_duplicates("sOTU").set_index("sOTU")
C["amphib_prevalence_pct"]=p90.map(sd.prevalence_pct)
M={}
for g in tgt:
    top=(C.blast_top_genus==g)&(C.blast_pct_identity>=97)&(C.aln_length>=100)
    tied=C.tied_genera.fillna("").str.split(";").apply(lambda x:g in x)&(C.blast_pct_identity>=97)
    ex=C.exact120_ref_genera.str.split(";").apply(lambda x:g in x)
    exr=ex&(C.exact120_n_refs>=2); exu=exr&(C.exact120_ref_genera==g)
    m90=C.metaTaxonomy90==g
    M[g]={k:v.fillna(False).astype(bool) for k,v in {"BLAST97_top":top,"BLAST97_tied":tied,"exact120":ex,"exact120_min2refs":exr,"exact120_unique_min2refs":exu,"metaTaxonomy90_97":m90,"any_method":top|tied|ex|m90}.items()}
inc=Tt>=1000; print("excluded samples:",Tt[~inc].to_dict())
cols={};summ=[]
for g,d in M.items():
    gk=g.replace(" (other)","_other")
    for mk,mask in d.items():
        r=RB.loc[mask[mask].index].sum() if mask.any() else pd.Series(0.0,index=RB.columns)
        cols[f"{gk}__{mk}"]=r.round(5); ri=r[inc]
        summ.append(dict(genus=g,method=mk,n_ASVs=int(mask.sum()),prevalence_pct=round(100*(ri>0).mean(),1),n_samples=int((ri>0).sum()),mean_rel_pct=round(ri.mean(),4),median_rel_pct=round(ri.median(),4),max_rel_pct=round(ri.max(),3)))
SB=pd.DataFrame(summ); W2=pd.DataFrame(cols)
print(SB.pivot(index="genus",columns="method",values="mean_rel_pct").to_string())
print(SB.pivot(index="genus",columns="method",values="prevalence_pct").to_string())
print(SB.pivot(index="genus",columns="method",values="n_ASVs").to_string())
