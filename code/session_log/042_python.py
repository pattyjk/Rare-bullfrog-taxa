# session cell 42 | 2026-10-15 21:02 | env: python | language: python
# files written: []


A0=A.set_index("ASV_ID")
ex0=A0[A0.mismatches==0]
def tgts(s): return ";".join(g for g in s.split(";") if g in tg)
C=B.copy()
C["exact120_ref_genera"]=ex0.matched_genera.reindex(C.index).fillna("")
C["exact120_targets"]=C.exact120_ref_genera.apply(tgts)
C["exact120_ref_species"]=ex0.top_ref_species.reindex(C.index).fillna("")
# amphibian link via 90-bp prefix
seq=pd.Series(q)
C["sequence"]=seq.reindex(C.index)
pref=C.sequence.str[:90]
amph_prev=(rel>0).mean()*100
inamph=pref.isin(cnt.index)
C["amphib_sOTU_match90"]=inamph
C["amphib_prevalence_pct"]=[round(amph_prev.get(p,np.nan),2) if p in cnt.index else np.nan for p in pref]
# note: cnt only holds target-family sOTUs; check against full sOTU list from taxonomy file for presence
C["in_meta_taxonomy_file"]=pref.isin(seqset)
sel=C[(C.top_hit_target_97!="")|((C.tied_targets!="")&(C.blast_pct_identity>=97))|(C.exact120_targets!="")]
print(len(sel))
# summary per genus
summ=[]
for g in tg:
    s1=(C.top_hit_target_97==g); s2=C.tied_targets.str.split(";").apply(lambda x:g in x)&(C.blast_pct_identity>=97); s3=C.exact120_targets.str.split(";").apply(lambda x:g in x)
    s3u=C.exact120_ref_genera==g
    summ.append(dict(genus=g,BLAST_top_hit_97=int(s1.sum()),BLAST_tied_best_97=int(s2.sum()),exact120_ref_match=int(s3.sum()),exact120_unique_to_genus=int(s3u.sum()),any_method=int((s1|s2|s3).sum())))
SM=pd.DataFrame(summ); print(SM.to_string())
print("ASVs with no BLAST hit:", len(q)-len(B))
