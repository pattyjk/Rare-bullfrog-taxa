# session cell 68 | 2026-10-15 21:58 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/antiSMASH_BGC_summary.xlsx","sha256":"095ab5c647365671516e307af2a61d6f187aa36bfd10b8341443dd60367c10ca"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/antiSMASH_regions.csv","sha256":"df5677922e59684464155a015cb47598650694cb382056999fe0d05c32f195f9"}]


ids={l.split("\t")[0]:json.loads(l.split("\t")[1])["id"] for l in open("as/jobs.tsv")}
AS["antiSMASH_result_url"]=AS.isolate.map({iso_map[k]:f"https://antismash.secondarymetabolites.org/upload/{v}/index.html" for k,v in ids.items()})
def simcls(x): return "no MIBiG hit" if pd.isna(x) else ("high (>=75%)" if x>=75 else ("medium (50-74%)" if x>=50 else "low (<50%)"))
AS["KCB_similarity_class"]=AS.KCB_similarity_pct.apply(simcls)
order_iso=["BF7","BF3","BF6","BF5","BF4","BF2","BF1"]
summ=[]
for iso in order_iso:
    s=AS[AS.isolate==iso]
    summ.append(dict(isolate=iso,taxonomy_MiGA=taxa[iso],n_regions=len(s),total_BGC_kb=round(s.length_kb.sum(),1),on_contig_edge=int(s.on_contig_edge.sum()),
        high_similarity_known=int((s.KCB_similarity_pct>=75).sum()),medium=int(s.KCB_similarity_pct.between(50,74.99).sum()),
        low_or_none_putatively_novel=int(((s.KCB_similarity_pct<50)|s.KCB_similarity_pct.isna()).sum()),
        region_types="; ".join(f"{k} ({v})" for k,v in s.products.str.split(";").explode().value_counts().items()),
        known_hits_ge50="; ".join(f"{r.KCB_top_hit} ({int(r.KCB_similarity_pct)}%)" for r in s[s.KCB_similarity_pct>=50].itertuples())))
SUM=pd.DataFrame(summ)
mat=AS.assign(p=AS.products.str.split(";")).explode("p").pivot_table(index="p",columns="isolate",values="region",aggfunc="count",fill_value=0)[order_iso]
notes3=pd.DataFrame({"note":[
"antiSMASH 7.1.0, public web server (antismash.secondarymetabolites.org), submitted 2026-10-16.",
"Input: contigs >=1000 bp (contig counts and genome sizes match the MiGA reports). Gene calling with Prodigal; detection strictness 'relaxed'.",
"Extra modules: KnownClusterBlast, MIBiG ClusterCompare, Cluster Pfam, ActiveSiteFinder, RREFinder, TTA, Pfam2GO. ClusterBlast and SubClusterBlast were not run.",
"KCB_similarity_pct is the share of genes in the top MIBiG reference cluster that have a significant BLAST hit in the region (as shown in the antiSMASH 'Similarity' column). It does not guarantee the same compound.",
"on_contig_edge = TRUE means the region sits at a contig end and may be incomplete. Fragmented assemblies (BF3: 222 contigs) split BGCs, so region counts can be inflated and similarities underestimated.",
"Hybrid regions (e.g. 'NRPS;NRP-metallophore') are counted once per product type in the type_matrix sheet.",
"Result pages are hosted on the antiSMASH server only for a limited time. The downloaded zips are saved as artifacts (open index.html)."]})
with pd.ExcelWriter("antiSMASH_BGC_summary.xlsx") as w:
    SUM.to_excel(w,sheet_name="per_isolate",index=False); AS.to_excel(w,sheet_name="all_regions",index=False)
    mat.to_excel(w,sheet_name="type_matrix"); notes3.to_excel(w,sheet_name="methods_notes",index=False)
AS.to_csv("antiSMASH_regions.csv",index=False)
print(SUM.drop(columns=["region_types"]).to_string()); print(mat.to_string())
