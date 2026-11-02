# session cell 51 | 2026-10-15 21:23 | env: python | language: python
# files written: []


g="Kineococcus"
seqs={};h=None
for l in open(f"meta/ref/{g}.fasta"):
    l=l.strip()
    if l.startswith(">"): h=l[1:]; seqs[h]=[]
    elif h: seqs[h].append(l.upper())
k90={};n120=0
for h,s in seqs.items():
    s="".join(s)
    for qq in (s,s.translate(rc)[::-1]):
        mm=pr.search(qq)
        if mm and len(qq)>=mm.end()+90:
            f90=qq[mm.end():mm.end()+90]
            if set(f90)<=set("ACGT"): k90[h]=f90; v2g.setdefault(f90,set()).add(g)
            if len(qq)>=mm.end()+120:
                fr=qq[mm.end():mm.end()+120]
                if set(fr)<=set("ACGT"):
                    sp=" ".join(h.split()[1:3]); ref120.setdefault(fr,{}).setdefault(g,collections.Counter())[sp]+=1; n120+=1
            break
refs[g]=k90
print("Kineococcus refs: 90bp",len(k90),"unique",len(set(k90.values())),"| 120bp",n120)
print("90bp variants shared with other genera:",sum(1 for v in set(k90.values()) if len(v2g[v])>1))
# bullfrog: BLAST
tg2=tg+["Kineococcus"]
best_k=best.groupby("qseqid").genus.apply(lambda x:"Kineococcus" in set(x))
Cx["blast_top_genus"]=B.blast_top_genus.reindex(ids_all)
blast97=((Cx.blast_top_genus=="Kineococcus")&(Cx.blast_pct_identity>=97)&(Cx.aln_length>=100)).fillna(False)
tied=(best_k.reindex(ids_all).fillna(False)&(Cx.blast_pct_identity>=97)).fillna(False).astype(bool)
# exact120
def ex_info(s):
    d=ref120.get(s,{}); return d
exk=pd.Series({a:("Kineococcus" in ex_info(q[a])) for a in ids_all})
nref=pd.Series({a:sum(sum(c.values()) for c in ex_info(q[a]).values()) for a in ids_all})
exr=exk&(nref>=2)
exu=exr&pd.Series({a:set(ex_info(q[a]))=={"Kineococcus"} for a in ids_all})
# meta taxonomy 90
ti["genus_tax"]=ti.Taxonomy.fillna("").str.rstrip("; ").str.split("; ").str[-1]
m90=p90.apply(lambda s: s in ti.index and ti.loc[s].genus_tax=="Kineococcus" and ti.loc[s].Per_ident>=97)
M["Kineococcus"]={"BLAST97_top":blast97,"BLAST97_tied":tied,"exact120":exk,"exact120_min2refs":exr,"exact120_unique_min2refs":exu,"metaTaxonomy90_97":m90,"any_method":blast97|tied|exk|m90}
for mk,mask in M["Kineococcus"].items(): print(mk,int(mask.sum()))
kids=M["Kineococcus"]["any_method"]; kids=kids[kids].index
print(B.reindex(kids)[["blast_top_hit","blast_pct_identity","tied_genera"]].to_string())
print("Kineococcus in tax file lineages:",(ti.genus_tax=="Kineococcus").sum())
