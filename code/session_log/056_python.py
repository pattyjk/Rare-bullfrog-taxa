# session cell 56 | 2026-10-15 21:28 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/MiGA_isolate_summary_7.csv","sha256":"da6fdc2db7681bbbeaa749b1b580d56faef144060455e1b3c70844d7aa7a4153"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/MiGA_isolate_summary_7.xlsx","sha256":"2f7ac407ce297d850895b9b23da2a6bbad273adbb8057960bd36809810ffe909"}]


def parse(k,s):
    s2=re.sub(r"\n\d+/\d+/\d+, [^\n]*\n[^\n]*MiGA Online\nPage \d+ of \d+\nhttps://[^\n]*\n","\n",s)
    flat=re.sub(r"\s+"," ",s2)
    r={"file":k,"MiGA_dataset":one(r"\n (\S[^\n]*typemat)\n",s2),"date":one(r"^(\d+/\d+/\d+)",s)}
    head=s2.split("Ribosomal and transfer")[0]
    for rk,lab in zip("pcofgs",["phylum","class","order","family","genus","species"]):
        r[f"MiGA_taxonomy_{lab}"]=one(rf"\n{rk}: ([^\n]+)",head)
    dist=one(r"Distances (.*?) Taxonomic classi",flat); r["closest_relatives"]=dist
    hits=re.findall(r"(?:were|was|and) ([A-Z][^()]*?GCF \d+ \d+(?: \([^)]*\))?) \(([\d.]+)% (ANI|AAI)\)",dist or "")
    if hits: r["closest_relative"]=hits[0][0]; r["closest_metric"]=hits[0][2]; r["closest_value_%"]=float(hits[0][1])
    r["classification"]=one(r"Taxonomic classi\S*cation (The dataset.*?\)\.) root",flat)
    c=r["classification"] or ""
    r["most_likely_rank"]=one(r"probably belongs to the (\w+) ",c)
    r["probable_taxon"]=one(r"probably belongs to the \w+ (.*?) \(p- ?value: ([\d.]+)\)",c)
    r["probable_taxon_p"]=num(one(r"probably belongs to the \w+ .*? \(p- ?value: ([\d.]+)\)",c))
    r["novelty"]=one(r"Taxonomic novelty (The dataset.*?) P-values:",flat)
    r["MiGA_novelty_call"]=novel(r["novelty"])
    r["novelty_pvals"]=one(r"P-values: (root.*?dataset [\d.]+)",flat)
    r["genus_pvalue"]=num(one(r"genus [A-Z]\w+ \(p-value ([\d.]+)",flat))
    m=re.search(r"Ribosomal RNA \(rRNAs\) MiGA found (.*?) Classify",flat); r["rRNAs"]=m.group(1) if m else None
    for lab,pat in [("16S_count",r"SSU: (\d+)"),("16S_max_len_bp",r"Max length: ([\d,]+)"),("tRNA_count",r"TRNA count: (\d+)"),("tRNA_AA",r"TRNA AA: (\d+)"),
                    ("completeness_%",r"Completeness: ([\d.]+) %"),("contamination_%",r"Contamination: ([\d.]+) %"),("quality",r"Quality: ([\d.]+)"),
                    ("predicted_proteins",r"Predicted proteins: ([\d,]+)"),("avg_protein_len_aa",r"Average length: ([\d.]+)"),("coding_density_%",r"Coding density: ([\d.]+)"),
                    ("contigs",r"Contigs: ([\d,]+)"),("N50_bp",r"N50: ([\d,]+) bp"),("genome_size_bp",r"Total length: ([\d,]+) bp"),("longest_contig_bp",r"Longest sequence: ([\d,]+)"),("GC_%",r"G\+C content: ([\d.]+)")]:
        r[lab]=num(one(pat,flat))
    r["quality_label"]=one(r"(\w+ quality genome)",flat)
    m=re.search(r"Fraction of proteins shared \(%\) (.*?) ([\d.]+) ([\d.]+) \(estimated\)",flat)
    if m: r["top_AAI_hit"]=m.group(1); r["top_AAI_%"]=float(m.group(2))
    return r
df7=pd.DataFrame([parse(k,s) for k,s in fullN.items()])
df7.loc[df7.file=="BF4","MiGA_dataset"]="BF4 typemat"
interp["BF3"]="Putative novel species in Kineosporiaceae, most likely Kineococcus (only AAI hits: 74.21% to K. rhizosphaerae, 74.04% to K. endophyticus; genus not significant, p=0.111). Folder label 'Kineo' matches."
df7["interpretation"]=df7.file.map(interp)
df7["rRNA_tRNA_note"]=df7.tRNA_count.apply(lambda x:"" if pd.notna(x) else "rRNA/tRNA section empty in report")
# check reproducibility vs earlier parse for the 6 repeats
for k in df7.file:
    if k in df.file.values:
        a=df7.set_index("file").loc[k]; b=df.set_index("file").loc[k]
        diffs=[c for c in ["closest_value_%","completeness_%","contamination_%","genome_size_bp","GC_%","contigs","N50_bp","predicted_proteins"] if a[c]!=b[c]]
        if diffs: print(k,"DIFF",diffs)
cols7=["file","MiGA_dataset","date","quality_label","MiGA_taxonomy_phylum","MiGA_taxonomy_class","MiGA_taxonomy_order","MiGA_taxonomy_family","MiGA_taxonomy_genus","MiGA_taxonomy_species",
"closest_relative","closest_metric","closest_value_%","top_AAI_hit","top_AAI_%","most_likely_rank","probable_taxon","probable_taxon_p","genus_pvalue","MiGA_novelty_call","interpretation",
"completeness_%","contamination_%","quality","genome_size_bp","contigs","N50_bp","longest_contig_bp","GC_%","predicted_proteins","avg_protein_len_aa","coding_density_%",
"16S_count","16S_max_len_bp","tRNA_count","tRNA_AA","rRNAs","rRNA_tRNA_note","closest_relatives","classification","novelty","novelty_pvals"]
o7=df7[cols7]
o7.to_excel("MiGA_isolate_summary_7.xlsx",index=False); o7.to_csv("MiGA_isolate_summary_7.csv",index=False)
print(o7[["file","MiGA_taxonomy_family","MiGA_taxonomy_genus","closest_relative","closest_metric","closest_value_%","MiGA_novelty_call","genome_size_bp","GC_%","contigs","completeness_%","contamination_%"]].to_string())
