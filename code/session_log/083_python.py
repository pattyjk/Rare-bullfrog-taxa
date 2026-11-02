# session cell 83 | 2026-10-17 17:50 | env: python | language: python
# files written: []


print("median percentile rank of matched present sOTUs:", a.percentile_rank.median(), "| top-100 sOTUs matched by AmphiBac:", int((present.mean_rel_pct.rank(ascending=False)<=100)[present.index.isin(agg.index)].sum()))
inh_s=agg.index[agg.function_class.str.startswith("Inhib")&agg.index.isin(RM.index)]
non_s=agg.index[(agg.function_class=="Non-inhibitory/mixed")&agg.index.isin(RM.index)]
fac_s=agg.index[agg.function_class.str.startswith("Facil")&agg.index.isin(RM.index)]
anyb=agg.index[agg.index.isin(RM.index)]
share=pd.DataFrame({"AmphiBac_matched_pct":RM.loc[anyb].sum(),"Bd_inhibitory_pct":RM.loc[inh_s].sum(),"non_inhibitory_pct":RM.loc[non_s].sum(),"facilitating_pct":RM.loc[fac_s].sum()})
share=share.join(mp[["host_order","host_species","lifestage","country","biogeo_realm"]])
print(share[["AmphiBac_matched_pct","Bd_inhibitory_pct","non_inhibitory_pct","facilitating_pct"]].describe().round(2).loc[["mean","50%","max"]].to_string())
print(share.groupby("host_order")[["AmphiBac_matched_pct","Bd_inhibitory_pct"]].median().round(2).to_string())
bf=share[share.host_species.str.contains("catesbeian",na=False)]; print("bullfrog samples in meta:",len(bf), bf[["AmphiBac_matched_pct","Bd_inhibitory_pct"]].median().round(2).to_dict())
tg7=["Rhodococcoides","Rhodococcus","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter","Sphingomonas","Kineococcus","Acinetobacter","Bacillus","Microbacterium"]
G=agg[agg.genus_meta.isin(tg7)]
print(G.groupby("genus_meta").agg(sOTUs=("n_isolates","size"),isolates=("n_isolates","sum"),inhib=("n_inhib","sum"),non=("n_non","sum"),max_prev=("prevalence_pct","max"),max_mean=("mean_rel_pct","max")).to_string())
