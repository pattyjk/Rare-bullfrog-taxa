# session cell 184 | 2026-11-02 11:13 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta_sOTU_stats.csv.gz","sha256":"d84e3d4d3f2abfd60ca346c1687e81e2be6306d6b67c00d6b68da5e3bf241bf8"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphibac_isolate_sOTU_candidates.csv.gz","sha256":"4d76d7e53335366bc08283baf56f7f7baccde21c93f0e38ee21cf063341b7d63"}]

out=I[['FastaID','Reference','Bd_function','inhib_score','has_V4b','n_exact','exact_sOTUs','nearest_mm','n_near','near_sOTUs','n_cand','set_orig']].rename(columns={'has_V4b':'has_V4'})
out.to_csv('amphibac_isolate_sOTU_candidates.csv.gz',index=False)
st.to_csv('meta_sOTU_stats.csv.gz',index=False)
b=I[I.Bd_function.notna()&(I.n_cand>0)]
print('Bd-tested matched isolates',len(b))
print('exact multi',(b.n_exact>1).sum(),' exact single',(b.n_exact==1).sum(),' near single',((b.n_exact==0)&(b.n_near==1)).sum(),' near tied',((b.n_exact==0)&(b.n_near>1)).sum())
nt=b[(b.n_exact==0)&(b.n_near>1)]
print('near-tie size median/max',nt.n_near.median(),nt.n_near.max())
# class change most-abundant vs rarest
ca=[clsd[most(c,max)] for c in nt.cand]; cr=[clsd[most(c,min)] for c in nt.cand]
print(pd.crosstab(pd.Series(ca,name='most_abund'),pd.Series(cr,name='rarest')))
em=b[b.n_exact>1]; print('exact-multi class mixes', Counter(tuple(sorted(set(clsd[s] for s in c.split(';')))) for c in em.cand).most_common(5))
# sOTU label conflict
Eo=E; g=Eo.groupby('sOTU').Bd_function
nl=g.nunique(); n=g.size()
print('sOTUs with >=2 isolates',(n>=2).sum(),' conflicting',((n>=2)&(nl>1)).sum(), round(100*((n>=2)&(nl>1)).sum()/(n>=2).sum(),1))
print('conflict by class', (pd.DataFrame({'n':n,'nl':nl,'cls':T.cls})[n>=2].assign(c=lambda d:d.nl>1).groupby('cls').c.agg(['size','mean'])).reindex(CL).round(2).to_dict())