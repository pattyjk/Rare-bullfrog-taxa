# session cell 6 | 2026-10-15 19:41 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/MiGA_isolate_summary.csv","sha256":"0393fc0da6411fecf767d744284dd9916545ce928f0939f87eacbbe6a9f06e7d"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/MiGA_isolate_summary.xlsx","sha256":"e7a909e8670c19af4d487a7c7a8eee1957059acf0d46e237c1c66f1b1d084a18"}]


df.loc[df.file=="BF4","MiGA_dataset"]="BF4 typemat"
def novel(t):
    if not t: return None
    m=re.search(r"(?:most likely|probably) belongs to an? (\w+) not represented in the database \(p- ?value: ([\d.]+)\)(?:(?!It probably).)*$",t) 
    hits=re.findall(r"belongs to an? (\w+) not represented in the database \(p- ?value: ([\d.]+)\)",t)
    # pick highest rank among hits
    order=["domain","phylum","class","order","family","genus","species","subspecies"]
    hits=sorted(hits,key=lambda h: order.index(h[0]) if h[0] in order else 99)
    return f"novel {hits[0][0]} (p={hits[0][1]})" if hits else None
df["MiGA_novelty_call"]=df.novelty.apply(novel)
df["rRNA_tRNA_note"]=df["tRNA_count"].apply(lambda x: "" if pd.notna(x) else "rRNA/tRNA section empty in report")
interp={
"HB58":"Acinetobacter oleivorans; borderline (ANI 95.25%, just above the 95% species cutoff; MiGA novelty test also suggests novel species p=0.019). Folder label 'Kineo' does not match.",
"BF4old":"Bacillus velezensis (ANI 98.81%). Folder label 'Caldifonti' does not match.",
"BF7":"Rhodococcoides fascians (syn. Rhodococcus fascians) (ANI 98.05%).",
"BF1":"Sphingomonas paucimobilis (ANI 99.49%).",
"BF2":"Putative novel species, possibly novel genus, in Xanthobacteraceae; only AAI hits (70.15% to Ancylobacter pratisalsi), genus not supported (p=0.122). Starkeya is in the same family.",
"BF246":"Putative novel Microbacterium species (ANI 85.94%, AAI 86.74% to M. mcarthurae).",
"BF4":"Putative novel Cryobacterium species (ANI 90.97% to C. soli; MiGA species p=0.046 conflicts, but ANI is below 95%).",
"BF5":"Putative novel Arthrobacter species (ANI 88.16% to A. humicola). CONFLICT: the other BF5 bin falls in Cryobacterium (different family).",
"BF6":"Putative novel species, genus Psychromicrobium (ANI 83.3%, AAI 84.71% to P. lacuslunae). Folder label 'Renibacterium' is the same family (Micrococcaceae).",
}
df["interpretation"]=df.file.map(interp)
cols=["file","MiGA_dataset","date","quality_label","MiGA_taxonomy_phylum","MiGA_taxonomy_class","MiGA_taxonomy_order","MiGA_taxonomy_family","MiGA_taxonomy_genus","MiGA_taxonomy_species",
"closest_relative","closest_metric","closest_value_%","top_AAI_hit","top_AAI_%","most_likely_rank","probable_taxon","probable_taxon_p","MiGA_novelty_call","interpretation",
"completeness_%","contamination_%","quality","genome_size_bp","contigs","N50_bp","longest_contig_bp","GC_%","predicted_proteins","avg_protein_len_aa","coding_density_%",
"16S_count","16S_max_len_bp","tRNA_count","tRNA_AA","rRNAs","rRNA_tRNA_note","closest_relatives","classification","novelty","novelty_pvals"]
out=df[cols]
out.to_csv("MiGA_isolate_summary.csv",index=False)
try:
    out.to_excel("MiGA_isolate_summary.xlsx",index=False); print("xlsx ok")
except Exception as e: print("no xlsx",e)
print(out[["file","MiGA_novelty_call"]].to_string())
