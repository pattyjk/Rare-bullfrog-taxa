# session cell 12 | 2026-10-15 19:48 | env: python | language: python
# files written: []

print(tax[tax.Hit_name.str.contains("kroppenstedtii",na=False)][["Hit_name","Per_ident","Taxonomy","hit_genus","hit_epithet","target"]].to_string())