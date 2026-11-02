# session cell 186 | 2026-11-02 11:14 | env: python | language: python
# files written: []

To,Eo=taxa('set_orig')
add('Majority of isolates inhibitory (used so far)','B. sOTU label rule (original tie rule)',To,To.fi>=0.5)
add('Any isolate inhibitory (AmphiBac subset-mode logic)','B. sOTU label rule (original tie rule)',To,To.fi>0)
add('All isolates inhibitory','B. sOTU label rule (original tie rule)',To,To.fi==1)
un=Eo.groupby('sOTU').Bd_function.nunique()==1
Tu=To[un.reindex(To.index).values]
add('Unanimous sOTUs only (conflicting sOTUs dropped)','B. sOTU label rule (original tie rule)',Tu,Tu.fi==1)
Tm=To[(Eo.groupby('sOTU').size()>=2).reindex(To.index).values]
add('Majority, sOTUs with >=2 isolates only','B. sOTU label rule (original tie rule)',Tm,Tm.fi>=0.5)
# top-hit analogue: one random isolate per sOTU
rr2=[]
grp=Eo.groupby('sOTU').Bd_function.apply(list)
for it in range(500):
    lab=pd.Series({s:(rng.choice(v)=='Inhibitory') for s,v in grp.items()}).reindex(To.index)
    rr2.append(summarise(To,lab))
rr2=pd.DataFrame(rr2)
d={k:rr2[k].median() for k in rr2}; d.update(scenario='One random isolate per sOTU (top-hit analogue; median of 500 draws)',group='B. sOTU label rule (original tie rule)'); res.append(d)
print(rr2[['inh_D','inh_I','inh_R','chi2_p']].quantile([.025,.5,.975]).round(3))
R_=pd.DataFrame(res)[['group','scenario','n_D','n_I','n_R','inh_D','inh_I','inh_R','chi2_p','spearman_abund_vs_inh']]
pd.set_option('display.width',250); print(R_.round(3).to_string(index=False))
# continuous score, each tie rule
sc=[]
for col,name in [('set_orig','most abundant'),('set_rare','rarest'),('set_all','all tied'),('set_unamb','drop tied'),('set_exact','exact only')]:
    Tt,_=taxa(col); Tt=Tt[Tt.med.notna()]
    sc.append(dict(tie_rule=name,median_D=Tt[Tt.cls=='Dominant'].med.median(),median_I=Tt[Tt.cls=='Intermediate'].med.median(),median_R=Tt[Tt.cls=='Rare'].med.median(),
                   kruskal_p=ss.kruskal(*[Tt[Tt.cls==c].med for c in CL])[1],spearman_abund_vs_score=ss.spearmanr(Tt.index.map(ab),Tt.med)[0],spearman_p=ss.spearmanr(Tt.index.map(ab),Tt.med)[1]))
SC=pd.DataFrame(sc); print(SC.round(3).to_string(index=False))