# session cell 67 | 2026-10-15 21:58 | env: python | language: python
# files written: []


iso_map={"BF1":"BF1","BF2":"BF2","BF3":"BF3","BF5":"BF5","BF4":"BF4","BF6":"BF6","BF7":"BF7"}
taxa={"BF1":"Sphingomonas paucimobilis","BF2":"Xanthobacteraceae sp. nov.","BF3":"Kineococcus sp. nov. (likely)","BF5":"Arthrobacter sp. nov.","BF4":"Cryobacterium sp. nov.","BF6":"Psychromicrobium sp. nov.","BF7":"Rhodococcoides fascians"}
rows=[]; vers={}
for d,iso in iso_map.items():
    jf=glob.glob(f"as/res/{d}/*.json")[0]; J=json.load(open(jf)); vers[iso]=J["version"]
    for ri,rec in enumerate(J["records"],1):
        if not rec.get("areas"): continue
        regf=[f for f in rec["features"] if f["type"]=="region"]
        kcb={r["region_number"]:r for r in rec["modules"].get("antismash.modules.clusterblast",{}).get("knowncluster",{}).get("results",[])}
        for ai,(area,rf) in enumerate(zip(rec["areas"],regf),1):
            q_=rf["qualifiers"]; rn=int(q_.get("region_number",[ai])[0])
            top=None
            k=kcb.get(rn)
            if k and k["ranking"]: top=k["ranking"][0]
            cats=sorted({p["category"] for p in area["protoclusters"].values()})
            rows.append(dict(isolate=iso,taxonomy_MiGA=taxa[iso],region=f"{ri}.{rn}",contig=rec["id"],start=area["start"]+1,end=area["end"],length_kb=round((area["end"]-area["start"])/1000,1),
                products=";".join(area["products"]),category=";".join(cats),on_contig_edge=q_.get("contig_edge",["?"])[0]=="True",
                KCB_top_hit=top[0]["description"] if top else "",KCB_accession=top[0]["accession"] if top else "",KCB_type=top[0]["cluster_type"] if top else "",
                KCB_similarity_pct=top[1]["similarity"] if top else np.nan))
AS=pd.DataFrame(rows)
print(vers); print(AS.groupby("isolate").size().to_dict())
print(AS[["isolate","region","products","length_kb","on_contig_edge","KCB_top_hit","KCB_similarity_pct"]].to_string())
