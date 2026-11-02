# session cell 111 | 2026-10-17 21:33 | env: python | language: python
# files written: []


mg=mg.copy()
U=mg.groupby(["sOTUs","hsp"]).agg(hg=("hg","first"),hf=("hf","first"),home_level=("home_level","first"),n_isolates=("FastaID","size"),
    n_inhib=("Bd_function",lambda x:(x=="Inhibitory").sum()),n_bd=("Bd_function",lambda x:x.notna().sum()),
    host_label=("Host.Species","first"),studies=("Reference",lambda x:"; ".join(sorted(set(map(str,x))))[:200])).reset_index()
U["function"]=np.where(U.n_bd==0,"untested",np.where(U.n_inhib/U.n_bd.clip(lower=1)>=0.5,"Inhibitory","Non-inhibitory/facilitating"))
hsp_arr=mpg.hsp.values; hg_arr=mpg.hg.values; hf_arr=mpg.hf.values
spN=mpg.hsp.value_counts(); okspp=spN[spN>=5].index
rows=[]
for _,u in U.iterrows():
    if u.home_level=="no relative in meta-analysis": continue
    ab=RM.loc[u.sOTUs.split(";")].sum().values
    home={"same host species":hsp_arr==u.hsp,"same host genus":hg_arr==u.hg,"same host family":hf_arr==u.hf}[u.home_level]
    oth=~home
    r=dict(home_n=int(home.sum()),home_mean=ab[home].mean(),home_prev=100*(ab[home]>0).mean(),other_mean=ab[oth].mean(),other_prev=100*(ab[oth]>0).mean(),overall_mean=ab.mean())
    # rank of home among host species (>=5 samples)
    sm=pd.Series(ab,index=mpg.index).groupby(mpg.hsp).mean().loc[okspp]
    if u.home_level=="same host species" and u.hsp in sm.index:
        r["home_rank_among_species"]=int(sm.rank(ascending=False,method="min")[u.hsp]); r["n_host_species"]=len(sm)
        r["home_is_top_host"]=r["home_rank_among_species"]==1
    r["top_host_species"]=sm.idxmax()
    rows.append(r)
R=pd.concat([U[U.home_level!="no relative in meta-analysis"].reset_index(drop=True),pd.DataFrame(rows)],axis=1)
pc=1e-4
R["log2_home_vs_other"]=np.log2((R.home_mean+pc)/(R.other_mean+pc))
print("units:",len(R), R.home_level.value_counts().to_dict())
for lv in ["same host species","same host genus","same host family"]:
    s=R[R.home_level==lv]; w=sst.wilcoxon(s.home_mean,s.other_mean)
    print(f"{lv}: n={len(s)} median home {s.home_mean.median():.4f}% vs other {s.other_mean.median():.4f}% | home>other {100*(s.home_mean>s.other_mean).mean():.0f}% | median log2 {s.log2_home_vs_other.median():.2f} | Wilcoxon p={w.pvalue:.2g} | prev {s.home_prev.median():.1f} vs {s.other_prev.median():.1f}")
sp=R[R.home_level=="same host species"]
print("home host is top host:",int(sp.home_is_top_host.sum()),"of",int(sp.home_is_top_host.notna().sum()), "median percentile rank:", round((100*sp.home_rank_among_species/sp.n_host_species).median(),1))
for fc in ["Inhibitory","Non-inhibitory/facilitating"]:
    s=R[R.function==fc]; print(fc, len(s), "median log2",round(s.log2_home_vs_other.median(),2), "home>other %",round(100*(s.home_mean>s.other_mean).mean(),1))
print("MWU log2 inhib vs non:",sst.mannwhitneyu(R[R.function=="Inhibitory"].log2_home_vs_other,R[R.function=="Non-inhibitory/facilitating"].log2_home_vs_other).pvalue)
