# session cell 185 | 2026-11-02 11:14 | env: python | language: python
# files written: []

rng=np.random.default_rng(1)
def summarise(T,lab):
    t=T.copy(); t['inh']=lab
    k=pd.crosstab(t.cls,t.inh).reindex(CL).fillna(0)
    p=ss.chi2_contingency(k)[1] if k.shape[1]==2 else np.nan
    pct=(100*t.groupby('cls').inh.mean()).reindex(CL)
    rho=ss.spearmanr(t.index.map(ab),t.inh.astype(float))[0]
    return dict(n_D=int((t.cls=='Dominant').sum()),n_I=int((t.cls=='Intermediate').sum()),n_R=int((t.cls=='Rare').sum()),
                inh_D=pct['Dominant'],inh_I=pct['Intermediate'],inh_R=pct['Rare'],chi2_p=p,spearman_abund_vs_inh=rho)
res=[]
def add(name,group,T,lab): d=summarise(T,lab); d.update(scenario=name,group=group); res.append(d)
# --- tie rules (majority label)
I['set_rare']=[c if ne>0 else (most(c,min) if c else '') for c,ne in zip(I.cand,I.n_exact)]
I['set_all']=I.cand
I['set_unamb']=np.where(I.n_cand==1,I.cand,'')
I['set_exact']=np.where(I.n_exact>0,I.cand,'')
for col,name in [('set_orig','Tied nearest -> most abundant sOTU (used so far)'),('set_rare','Tied nearest -> rarest sOTU'),
                 ('set_all','Tied nearest -> all tied sOTUs'),('set_unamb','Drop isolates with tied nearest sOTUs'),('set_exact','Exact matches only (0 mismatches)')]:
    Tt,_=taxa(col); add(name,'A. isolate-to-sOTU tie rule (majority label)',Tt,Tt.fi>=0.5)
# random tie draws
rr=[]
for it in range(500):
    I['set_rand']=[c if ne>0 else (rng.choice(c.split(';')) if c else '') for c,ne in zip(I.cand,I.n_exact)]
    Tt,_=taxa('set_rand'); rr.append(summarise(Tt,Tt.fi>=0.5))
rr=pd.DataFrame(rr)
d={k:rr[k].median() for k in rr}; d.update(scenario='Tied nearest -> random sOTU (median of 500 draws)',group='A. isolate-to-sOTU tie rule (majority label)'); res.append(d)
print(rr[['inh_D','inh_I','inh_R','chi2_p']].quantile([.025,.5,.975]).round(3))