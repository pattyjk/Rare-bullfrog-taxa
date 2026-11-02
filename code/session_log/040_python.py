# session cell 40 | 2026-10-15 21:01 | env: python | language: python
# files written: []


b=pd.read_csv("meta/asv_blast.tsv",sep="\t",header=None,names="qseqid sacc pident length mismatch gapopen qstart qend evalue bitscore staxids sscinames stitle".split())
print("Rhodococcoides names in db:", b.sscinames.str.contains("Rhodococcoides").sum(), "Psychromicrobium:", b.sscinames.str.contains("Psychromicrobium").sum())
def gen_of(name):
    w=str(name).split()
    if len(w)>=2 and w[0] in("Rhodococcus","Rhodococcoides") and w[1] in rhodoides: return "Rhodococcoides"
    if w and w[0]=="Rhodococcus": return "Rhodococcus (other)"
    return w[0] if w else None
b["genus"]=b.sscinames.apply(gen_of)
best=b[b.bitscore==b.groupby("qseqid").bitscore.transform("max")]
first=b.groupby("qseqid",sort=False).head(1).set_index("qseqid")
agg=best.groupby("qseqid").agg(tied_genera=("genus",lambda x:";".join(sorted(set(x)))),n_tied_hits=("sacc","size"),
                               tied_species=("sscinames",lambda x:"; ".join(pd.Series(x).value_counts().index[:5])))
B=first[["sscinames","genus","pident","length"]].rename(columns={"sscinames":"blast_top_hit","genus":"blast_top_genus","pident":"blast_pct_identity","length":"aln_length"}).join(agg)
tg=tgt6+["Rhodococcus (other)"]
B["top_hit_target_97"]=np.where((B.blast_pct_identity>=97)&(B.aln_length>=100)&B.blast_top_genus.isin(tg),B.blast_top_genus,"")
B["tied_targets"]=B.tied_genera.apply(lambda s:";".join(g for g in s.split(";") if g in tg))
print(B.top_hit_target_97.value_counts().to_dict())
print((B[B.blast_pct_identity>=97].tied_targets!="").sum())
print(B[B.tied_targets.str.contains("Rhodococcoides|Psychromicrobium|Cryobacterium")][["blast_top_hit","blast_pct_identity","tied_genera","n_tied_hits"]].to_string())
