# session cell 4 | 2026-10-15 19:40 | env: python | language: python
# files written: []


import re, pandas as pd
def one(pat, s, g=1, flags=re.S):
    m=re.search(pat,s,flags); return m.group(g).strip() if m else None
def num(x):
    return float(x.replace(",","")) if x else None
rows=[]
for k,s in full.items():
    s2=re.sub(r"\n\d+/\d+/\d+, [^\n]*\n[^\n]*MiGA Online\nPage \d+ of \d+\nhttps://[^\n]*\n","\n",s)  # drop page headers
    flat=re.sub(r"\s+"," ",s2)
    r={"file":k,"MiGA_dataset":one(r"\n (\S[^\n]*typemat)\n",s2),"date":one(r"^(\d+/\d+/\d+)",s)}
    tax={}
    for rk,lab in zip("dpcofgs",["domain","phylum","class","order","family","genus","species"]):
        tax[lab]=one(rf"\n{rk}: ([^\n]+)",s2.split("Ribosomal and transfer")[0])
    r.update({f"MiGA_taxonomy_{x}":tax[x] for x in ["phylum","class","order","family","genus","species"]})
    dist=one(r"Distances (.*?) Taxonomic classi",flat)
    r["closest_relatives"]=dist
    hits=re.findall(r"(?:were|was|and) ([A-Z][^()]*?GCF \d+ \d+(?: \([^)]*\))?) \(([\d.]+)% (ANI|AAI)\)",dist or "")
    if hits:
        r["closest_relative"]=hits[0][0]; r["closest_metric"]=hits[0][2]; r["closest_value_%"]=float(hits[0][1])
    r["classification"]=one(r"Taxonomic classi\S*cation (The dataset.*?\)\.) root",flat)
    r["most_likely_rank"]=one(r"probably belongs to the (\w+) ",r["classification"] or "")
    r["probable_taxon"]=one(r"probably belongs to the \w+ (.*?) \(p- ?value: ([\d.]+)\)",r["classification"] or "")
    r["probable_taxon_p"]=num(one(r"probably belongs to the \w+ .*? \(p- ?value: ([\d.]+)\)",r["classification"] or ""))
    r["novelty"]=one(r"Taxonomic novelty (The dataset.*?) P-values:",flat)
    r["novel_rank_p<=0.1"]=one(r"probably belongs to a (\w+) not represented",r["novelty"] or "")
    r["novelty_pvals"]=one(r"P-values: (root.*?dataset [\d.]+)",flat)
    m=re.search(r"Ribosomal RNA \(rRNAs\) MiGA found (.*?) Classify",flat); r["rRNAs"]=m.group(1) if m else None
    for lab,pat in [("16S_count",r"SSU: (\d+)"),("16S_max_len_bp",r"Max length: ([\d,]+)"),("tRNA_count",r"TRNA count: (\d+)"),("tRNA_AA",r"TRNA AA: (\d+)"),
                    ("completeness_%",r"Completeness: ([\d.]+) %"),("contamination_%",r"Contamination: ([\d.]+) %"),("quality",r"Quality: ([\d.]+)"),
                    ("predicted_proteins",r"Predicted proteins: ([\d,]+)"),("avg_protein_len_aa",r"Average length: ([\d.]+)"),("coding_density_%",r"Coding density: ([\d.]+)"),
                    ("contigs",r"Contigs: ([\d,]+)"),("N50_bp",r"N50: ([\d,]+) bp"),("genome_size_bp",r"Total length: ([\d,]+) bp"),("longest_contig_bp",r"Longest sequence: ([\d,]+)"),("GC_%",r"G\+C content: ([\d.]+)")]:
        r[lab]=num(one(pat,flat))
    r["quality_label"]=one(r"(\w+ quality genome)",flat)
    # top AAI table hit
    m=re.search(r"Fraction of proteins shared \(%\) (.*?) ([\d.]+) ([\d.]+) \(estimated\)",flat)
    if m: r["top_AAI_hit"]=m.group(1); r["top_AAI_%"]=float(m.group(2))
    rows.append(r)
df=pd.DataFrame(rows)
pd.set_option("display.width",250,"display.max_columns",50,"display.max_colwidth",60)
print(df[["file","MiGA_dataset","MiGA_taxonomy_genus","MiGA_taxonomy_species","closest_relative","closest_metric","closest_value_%","top_AAI_hit","top_AAI_%"]].to_string())
print(df[["file","most_likely_rank","probable_taxon","probable_taxon_p","novel_rank_p<=0.1","completeness_%","contamination_%","quality","contigs","N50_bp","genome_size_bp","GC_%","predicted_proteins","tRNA_count","16S_count"]].to_string())
