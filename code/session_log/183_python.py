# session cell 183 | 2026-11-02 11:13 | env: python | language: python
# files written: []

import scipy.stats as ss
for idx in amb.index:
    w=amb.at[idx,'w']; d=(M!=np.frombuffer(w.encode(),dtype=np.uint8)).sum(1); dm=int(d.min())
    I.at[idx,'nearest_mm']=dm
    if dm<=3: near=list(st.sOTU.values[d==dm]); I.at[idx,'near_sOTUs']=';'.join(near); I.at[idx,'n_near']=len(near)
def norm_bd(x):
    x=str(x).strip()
    if x in ("Inhibitory","inhibitory","Inhibitor"): return "Inhibitory"
    if x in ("NonInhibitory","NoEffect","No Effect","No effect","Non-inhibitor","ns"): return "Non-inhibitory"
    if x in ("Facilitating","enhancing"): return "Facilitating"
    return None
I=I.merge(meta[['FastaID','Bd.Inhibition','Proportional.Growth.Bd','Reference']],on='FastaID',how='left')
I['Bd_function']=I['Bd.Inhibition'].map(norm_bd)
I['inhib_score']=(1-pd.to_numeric(I['Proportional.Growth.Bd'],errors='coerce'))*100
# candidate set = all exact sOTUs, else all tied nearest
I['cand']=np.where(I.n_exact>0,I.exact_sOTUs,I.near_sOTUs)
I['n_cand']=np.where(I.n_exact>0,I.n_exact,I.n_near)
print('detected',(I.n_cand>0).sum())
def most(c,f): 
    L=c.split(';'); return f(L,key=lambda s:ab[s])
# original rule: exact -> all; near tie -> most abundant
I['set_orig']=[c if ne>0 else (most(c,max) if c else '') for c,ne in zip(I.cand,I.n_exact)]
CL=['Dominant','Intermediate','Rare']
cls=lambda a:'Dominant' if a>=1 else ('Intermediate' if a>=0.1 else 'Rare')
st['cls']=st.mean_rel_pct.map(cls); clsd=dict(zip(st.sOTU,st.cls))
def taxa(col,df=None,rule='majority'):
    df=I if df is None else df
    b=df[df.Bd_function.notna()&(df[col]!='')]
    E=b.assign(sOTU=b[col].str.split(';')).explode('sOTU')
    g=E.groupby('sOTU')
    T=pd.DataFrame({'n':g.FastaID.nunique(),'fi':g.Bd_function.apply(lambda x:(x=='Inhibitory').mean()),
                    'ff':g.Bd_function.apply(lambda x:(x=='Facilitating').mean()),'med':g.inhib_score.median()})
    T['cls']=T.index.map(clsd); return T,E
T,E=taxa('set_orig')
T['inh']=T.fi>=0.5
print(T.cls.value_counts().reindex(CL).tolist(), (100*T.groupby('cls').inh.mean()).reindex(CL).round(0).tolist())
print('isolates',E.FastaID.nunique())