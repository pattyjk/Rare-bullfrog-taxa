# session cell 11 | 2026-10-15 19:48 | env: python | language: python
# files written: []


for pat in ["fascians","corynebacterioides","kroppenstedtii","kyoton","navarretei","trifolii","yunnanens","Psychromicrob","lacuslunae","Rhodococcoides"]:
    print(pat, tax.Hit_name.str.contains(pat,na=False).sum())
print(tax[tax.target=="Arthrobacter"].Hit_name.str.split().str[:2].str.join(" ").value_counts().head(10).to_dict())
print(tax[tax.target=="Ancylobacter"].Hit_name.str.split().str[:2].str.join(" ").value_counts().head(5).to_dict())
# Micrococcaceae genera present
mic=tax[tax.Taxonomy.str.contains("Micrococcaceae",na=False)]
print(mic.genus_tax.value_counts().head(15).to_dict())
dup=tax[tax.sOTU.duplicated(keep=False)]; print("dup sOTUs w/ target:",dup[dup.target.notna()].sOTU.nunique())
