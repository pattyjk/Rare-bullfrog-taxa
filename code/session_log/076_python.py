# session cell 76 | 2026-10-17 16:58 | env: python | language: python
# files written: []


anyf=pd.DataFrame({g:d["any_method"] for g,d in M.items()})
C["target_genera_any_method"]=anyf.apply(lambda r:";".join(g for g,v in r.items() if v),axis=1)
C["methods"]=[ "; ".join(f"{g}: "+",".join(mk for mk,m in M[g].items() if m[a] and mk!="any_method") for g in tgt if M[g]["any_method"][a]) for a in C.index]
C["mean_rel_pct"]=RB.mean(axis=1).round(4); C["prevalence_pct"]=((RB>0).mean(axis=1)*100).round(1); C["max_rel_pct"]=RB.max(axis=1).round(3)
D=C[C.target_genera_any_method!=""].copy()
pd.set_option("display.width",250,"display.max_colwidth",70)
print(D[D.target_genera_any_method!="Sphingomonas"][["methods","blast_top_hit","blast_pct_identity","tied_genera","exact120_n_refs","exact120_ref_species","mean_rel_pct","prevalence_pct"]].to_string())
