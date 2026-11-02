# session cell 45 | 2026-10-15 21:02 | env: python | language: python
# files written: []


amph_prev=(rel>0).mean(axis=1)*100
C["amphib_prevalence_pct"]=[round(amph_prev[p],2) if p in amph_prev.index else np.nan for p in pref]
print("Rhodococcoides ASV amphib prevalence:", C.loc["5cb28bd0b2eea6847a6328712217c141","amphib_prevalence_pct"])
sup=[]
for sid in M2["Cryobacterium"]:
    nC=sum(1 for h,v in refs["Cryobacterium"].items() if v==sid)
    sp=collections.Counter(" ".join(h.split()[1:3]) for h,v in refs["Cryobacterium"].items() if v==sid)
    sup.append(dict(n_cryo_refs=nC,species=dict(sp.most_common(2)),prev=round(amph_prev[sid],2),blast=" ".join(str(ti.loc[sid].Hit_name).split()[:2]),ident=ti.loc[sid].Per_ident))
print(pd.DataFrame(sup).sort_values("prev",ascending=False).to_string())
