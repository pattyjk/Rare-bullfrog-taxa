# session cell 43 | 2026-10-15 21:02 | env: python | language: python
# files written: []


def tgts(s): return ";".join(("Rhodococcus (other)" if g=="Rhodococcus" else g) for g in s.split(";") if g in tg or g=="Rhodococcus")
C["exact120_targets"]=C.exact120_ref_genera.apply(tgts)
summ=[];flags={}
for g in tg:
    s1=(C.top_hit_target_97==g); s2=C.tied_targets.str.split(";").apply(lambda x:g in x)&(C.blast_pct_identity>=97); s3=C.exact120_targets.str.split(";").apply(lambda x:g in x)
    s3u=C.exact120_targets==g if g!="Rhodococcus (other)" else C.exact120_ref_genera=="Rhodococcus"
    flags[g]=s1|s2|s3
    summ.append(dict(genus=g,BLAST_top_hit_97=int(s1.sum()),BLAST_tied_best_97=int(s2.sum()),exact120_ref_match=int(s3.sum()),exact120_unique_among_10_genera=int(s3u.sum()),ASVs_any_method=int((s1|s2|s3).sum())))
SM=pd.DataFrame(summ); print(SM.to_string())
anyf=pd.DataFrame(flags).any(axis=1)
C["target_genera_any_method"]=pd.DataFrame(flags).apply(lambda r:";".join(g for g,v in r.items() if v),axis=1)
D2=C[anyf].copy()
cols2=["target_genera_any_method","top_hit_target_97","tied_targets","exact120_targets","blast_top_hit","blast_pct_identity","aln_length","tied_genera","n_tied_hits","tied_species","exact120_ref_genera","exact120_ref_species","in_meta_taxonomy_file","amphib_sOTU_match90","amphib_prevalence_pct","sequence"]
D2=D2[cols2].sort_values(["target_genera_any_method","blast_pct_identity"],ascending=[True,False]); D2.index.name="ASV_ID"
print(D2[~D2.target_genera_any_method.isin(["Sphingomonas","Ancylobacter"])][["target_genera_any_method","blast_top_hit","blast_pct_identity","tied_genera","exact120_ref_species","amphib_prevalence_pct"]].to_string())
print("in amphib:",D2.amphib_sOTU_match90.sum(),"of",len(D2))
