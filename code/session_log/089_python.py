# session cell 89 | 2026-10-17 17:52 | env: python | language: python
# files written: []


ab_=agg[agg.mean_rel_pct==0].sort_values("n_isolates",ascending=False)
print(ab_[["n_isolates","top_hit_meta","countries"]].head(8).to_string())
print("absent sOTUs in amphib table index at all:", ab_.index.isin(ST.index).sum(), "of", len(ab_))
s0=ab_.index[0]; print(s0[:30], s0 in ST.index, ST.loc[s0].to_dict() if s0 in ST.index else None)
