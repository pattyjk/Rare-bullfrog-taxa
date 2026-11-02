# session cell 187 | 2026-11-02 11:15 | env: python | language: python
# files written: []

print('isolates per sOTU (median, mean) by class:', To.groupby('cls').n.agg(['median','mean']).reindex(CL).round(1).to_dict())
iso=[];cul=[]
tot_cls=st.cls.value_counts().reindex(CL)
for col,name in [('set_orig','most abundant (used so far)'),('set_rare','rarest'),('set_all','all tied'),('set_unamb','drop tied'),('set_exact','exact only')]:
    b=I[I.Bd_function.notna()&(I[col]!='')].copy()
    if col=='set_all':  # isolate takes the class of its most abundant AND rarest? use fractional: each tied sOTU weight 1/k
        E2=b.assign(sOTU=b[col].str.split(';')).explode('sOTU'); E2['w']=1/E2.groupby('FastaID').sOTU.transform('size')
    else:
        E2=b.assign(sOTU=b[col]); E2['w']=1.0
    E2['cls']=E2.sOTU.map(clsd); E2['inh']=(E2.Bd_function=='Inhibitory')
    g=E2.groupby('cls').apply(lambda d:pd.Series({'n':d.w.sum(),'pct':100*(d.w*d.inh).sum()/d.w.sum()})).reindex(CL)
    iso.append(dict(tie_rule=name,n_D=g.n['Dominant'],n_I=g.n['Intermediate'],n_R=g.n['Rare'],inh_D=g.pct['Dominant'],inh_I=g.pct['Intermediate'],inh_R=g.pct['Rare']))
    a=I[I[col]!='']; cs=set(';'.join(a[col]).split(';'))
    c=st[st.sOTU.isin(cs)].cls.value_counts().reindex(CL).fillna(0)
    cul.append(dict(tie_rule=name,cultured_D=f"{int(c['Dominant'])}/{tot_cls['Dominant']}",cultured_I=f"{int(c['Intermediate'])}/{tot_cls['Intermediate']}",cultured_R=f"{int(c['Rare'])}/{tot_cls['Rare']}"))
ISO=pd.DataFrame(iso); CUL=pd.DataFrame(cul)
print(ISO.round(1).to_string(index=False)); print(CUL.to_string(index=False))