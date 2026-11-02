# session cell 112 | 2026-10-17 21:34 | env: python | language: python
# files written: []


rows=[]
for i,u in R.iterrows():
    ab=pd.Series(RM.loc[u.sOTUs.split(";")].sum().values,index=mpg.index)
    lvlcol={"same host species":"hsp","same host genus":"hg","same host family":"hf"}[u.home_level]
    key={"hsp":u.hsp,"hg":u.hg,"hf":u.hf}[lvlcol]
    grp=ab.groupby(mpg[lvlcol]); gm=grp.mean(); gn=grp.size(); gm=gm[gn>=5]
    if key not in gm.index: rows.append({}); continue
    others=gm.drop(key)
    rows.append(dict(home_group_mean=gm[key],median_other_groups=others.median(),pct_other_groups_below_home=100*(others<gm[key]).mean(),
                     n_groups=len(gm),home_group_rank=int(gm.rank(ascending=False,method="min")[key])))
H=pd.DataFrame(rows,index=R.index); R2=pd.concat([R,H],axis=1)
R2["log2_home_vs_median_other_group"]=np.log2((R2.home_group_mean+pc)/(R2.median_other_groups+pc))
for lv in ["same host species","same host genus","same host family"]:
    s=R2[(R2.home_level==lv)&R2.home_group_mean.notna()]
    w=sst.wilcoxon(s.home_group_mean,s.median_other_groups)
    print(f"{lv}: n={len(s)} | home > median other {100*(s.home_group_mean>s.median_other_groups).mean():.0f}% | median % of other groups below home {s.pct_other_groups_below_home.median():.0f}% | median log2 {s.log2_home_vs_median_other_group.median():.2f} | Wilcoxon p={w.pvalue:.2g} | top group: {int((s.home_group_rank==1).sum())}")
s=R2[R2.home_group_mean.notna()]
for fc in ["Inhibitory","Non-inhibitory/facilitating"]:
    t=s[s.function==fc]; print(fc,len(t),"median % other groups below home",round(t.pct_other_groups_below_home.median(),1),"home>median",round(100*(t.home_group_mean>t.median_other_groups).mean(),1))
print("MWU pct_below inhib vs non:",sst.mannwhitneyu(s[s.function=="Inhibitory"].pct_other_groups_below_home,s[s.function=="Non-inhibitory/facilitating"].pct_other_groups_below_home).pvalue)
# how many other groups have zero (ties at 0)?
print("median fraction of other host species where taxon absent:", round(s.apply(lambda r: np.nan,axis=1).isna().mean(),2))
