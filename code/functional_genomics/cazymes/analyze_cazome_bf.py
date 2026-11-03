"""Post-processing of run_cazome_bf.py output for BF1-BF7.
(1) full CAZome  (2) mucin O-glycan panel + catabolic-operon neighbourhoods  (3) fungal cell-wall panel + DeepSig
Thresholds (stated):
 - dbCAN: i-Evalue < 1e-15 if HMM alignment > 80 aa else < 1e-5; HMM coverage > 0.35 (applied in run_cazome_bf.py).
 - Overlap resolution within a protein: two domains overlapping by > 50% of the shorter envelope -> keep best i-Evalue.
 - Subfamilies (GH5_x, GH13_x, GT2_* ...) collapsed to family. Counts = distinct proteins with >=1 domain of the family.
 - Pfam GA, NCBIfam/TIGRFAM TC cutoffs.
 - Operon neighbourhood: genes on the same contig joined by single linkage when the intergenic gap
   (end of one CDS to start of the next) is <= 10,000 bp.
 - Signal peptide: DeepSig DNN class 'SignalPeptide' (gramn threshold 0.7575 for BF1/BF2, gramp 0.72 for BF3-BF7);
   Tat motif: TIGR01409 or PF10518 hit starting within residue 1-60.
 - Most-variable families: rank by log2((max+1)/(min+1)) of per-isolate protein counts, ties by range.
"""
import re, glob, numpy as np, pandas as pd, pyhmmer, os
from pyhmmer.plan7 import HMMFile
from pyhmmer.easel import SequenceFile, Alphabet
from substrate_map import family, cazy_class, substrate
ISO = [f"BF{i}" for i in range(1, 8)]
TAXON = {"BF1":"Sphingomonas paucimobilis","BF2":"Xanthobacteraceae sp. nov.","BF3":"Kineococcus sp. nov.","BF4":"Cryobacterium sp. nov.",
         "BF5":"Arthrobacter sp. nov.","BF6":"Psychromicrobium sp. nov.","BF7":"Rhodococcoides fascians"}
NPROT = {i: sum(1 for l in open(f"fg_inputs/{i}/{i}.faa") if l.startswith(">")) for i in ISO}
os.makedirs("out", exist_ok=True)
def wide(df, idx, col="n"):
    m = df.pivot_table(index=idx, columns="isolate", values=col, aggfunc="sum", fill_value=0)
    return m.reindex(columns=ISO, fill_value=0)

# ---------- GFF ----------
g = []
for i in ISO:
    for l in open(f"fg_inputs/{i}/{i}.gff"):
        if l.startswith("##FASTA"): break
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        g.append((i, a["ID"], f[0], int(f[3]), int(f[4]), f[6], a.get("gene", ""), re.sub("%2C", ",", a.get("product", ""))))
gff = pd.DataFrame(g, columns=["isolate","protein","contig","start","end","strand","prokka_gene","prokka_product"])
assert all(gff.groupby("isolate").size()[i] == NPROT[i] for i in ISO)

# ---------- (1) full CAZome ----------
raw = pd.read_csv("raw_dbcan_hits.tsv", sep="\t")
raw["family"] = raw.model.map(family); raw["class"] = raw.family.map(cazy_class)
keep = []
for (iso, p), d in raw.sort_values("i_evalue").groupby(["isolate","protein"], sort=False):
    kept = []
    for r in d.itertuples():
        ok = True
        for k in kept:
            ov = min(r.env_to, k.env_to) - max(r.env_from, k.env_from) + 1
            if ov > 0.5 * min(r.env_to - r.env_from + 1, k.env_to - k.env_from + 1): ok = False; break
        if ok: kept.append(r)
    keep += [k.Index for k in kept]
dom = raw.loc[keep].sort_values(["isolate","protein","env_from"]).copy()
dom = dom.merge(gff, on=["isolate","protein"], how="left")
dom.insert(1, "taxon", dom.isolate.map(TAXON))
dom["substrate"] = dom.family.map(substrate)
dom.to_csv("out/cazome_domain_hits.tsv", sep="\t", index=False)
print("raw domains", len(raw), "after overlap", len(dom), "non-CAZy modules:", sorted(set(dom.family[dom["class"]=="Other module"])))

fam = dom[dom["class"] != "Other module"].drop_duplicates(["isolate","protein","family"])
fam_long = fam.groupby(["isolate","class","family"]).size().rename("n").reset_index()
fam_long["per_1000_proteins"] = (fam_long.n / fam_long.isolate.map(NPROT) * 1000).round(3)
fam_long.to_csv("out/cazy_family_by_isolate_long.tsv", sep="\t", index=False)
fw = wide(fam_long, ["class","family"])
order = {c: k for k, c in enumerate(["GH","GT","PL","CE","AA","CBM"])}
fw = fw.reset_index(); fw["o"] = fw["class"].map(order); fw["num"] = fw.family.str.extract(r"(\d+)").astype(int)
fw = fw.sort_values(["o","num"]).drop(columns=["o","num"]).set_index(["class","family"])
fw.to_csv("out/cazy_family_by_isolate_matrix.tsv", sep="\t")

deg = fam[fam["class"].isin(["GH","PL","CE","AA","CBM"])]
sub = deg.drop_duplicates(["isolate","protein","substrate"]).groupby(["isolate","substrate"]).size().rename("n").reset_index()
sub.to_csv("out/substrate_by_isolate_long.tsv", sep="\t", index=False)
wide(sub, "substrate").to_csv("out/substrate_by_isolate_matrix.tsv", sep="\t")
gt = fam[fam["class"] == "GT"].groupby(["isolate","family"]).size().rename("n").reset_index()
wide(gt, "family").to_csv("out/glycosyltransferase_by_isolate_matrix.tsv", sep="\t")

tot = pd.DataFrame({"isolate": ISO, "taxon": [TAXON[i] for i in ISO], "proteins": [NPROT[i] for i in ISO]})
tot["degradative_CAZyme_proteins"] = [deg[deg.isolate==i].protein.nunique() for i in ISO]
tot["degradative_per_1000_proteins"] = (tot.degradative_CAZyme_proteins / tot.proteins * 1000).round(2)
for c in ["GH","GT","PL","CE","AA","CBM"]:
    tot[f"{c}_proteins"] = [fam[(fam.isolate==i)&(fam["class"]==c)].protein.nunique() for i in ISO]
tot["GT_per_1000_proteins"] = (tot.GT_proteins / tot.proteins * 1000).round(2)
tot["any_CAZyme_proteins"] = [fam[fam.isolate==i].protein.nunique() for i in ISO]
tot["distinct_degradative_families"] = [deg[deg.isolate==i].family.nunique() for i in ISO]
tot.to_csv("out/cazome_totals_by_isolate.tsv", sep="\t", index=False)
print(tot.to_string())

m = fw.copy(); mx, mn = m.max(1), m.min(1)
var = pd.DataFrame({"min": mn, "max": mx, "range": mx - mn, "log2_max1_min1": np.log2((mx + 1) / (mn + 1)).round(3),
                    "isolates_present": (m > 0).sum(1)}).join(m)
var = var.sort_values(["log2_max1_min1","range"], ascending=False)
var.to_csv("out/cazy_family_variability_ranked.tsv", sep="\t")
print(var.head(12).to_string())

# ---------- panel hits (Pfam/NCBIfam) ----------
LAB = {l.split("\t")[0].split(".")[0]: None for l in open("db/panel_hmm_versions.tsv")}
from importlib import util
PANEL = {"PF00884":"Sulfatase","PF03781":"FGE-sulfatase","PF13402":"Peptidase_M60","PF10462":"Peptidase_M66","PF02128":"Peptidase_M36","PF00026":"Asp","PF03629":"SASA","PF04131":"nanE",
 "TIGR00683":"nanA","NF003164":"nanA","NF003461":"nanK","NF047821":"nanK","NF002231":"nanE",
 "TIGR01089":"fucI","NF008220":"fucI","PF02952":"fucI_domain","PF07881":"fucI_domain","PF07882":"fucI_domain",
 "TIGR02628":"fucK","TIGR01086":"fucA","NF004527":"fucA","NF005302":"fucA","NF005984":"fucA","NF040649":"fucA",
 "TIGR02638":"fucO","NF007911":"fucO","NF011949":"fucU","PF05025":"RbsD_FucU(generic)","TIGR00885":"fucP","NF007524":"fucP",
 "TIGR00221":"nagA","NF008371":"nagA","TIGR00502":"nagB","NF001682":"nagB","NF001684":"nagB","NF001685":"nagB","NF009022":"nagB","NF041128":"nagB","NF046059":"nagB-II","PF01182":"Glucosamine_iso(generic)",
 "NF007291":"agaS","TIGR02815":"agaS","NF007288":"aga_PTS","NF007289":"aga_PTS","NF007359":"aga_PTS","NF008508":"aga_PTS","NF040757":"aga_PTS","NF040760":"aga_PTS","NF040761":"aga_PTS",
 "TIGR01858":"kbaY_gatY","NF009375":"kbaY_gatY","TIGR02810":"kbaZ_gatZ","NF012002":"kbaZ_gatZ","PF08013":"kbaZ_gatZ","TIGR01231":"lacC","NF010033":"lacC","TIGR01232":"lacD",
 "PF03480":"TRAP_DctP","PF06808":"TRAP_DctM","PF04290":"TRAP_DctQ","PF07690":"MFS_1","PF00083":"Sugar_tr",
 "PF01547":"ABC_SBP_1","PF13416":"ABC_SBP_8","PF00528":"ABC_permease_BPD","TIGR01409":"Tat_signal","PF10518":"Tat_signal"}
pan = pd.read_csv("raw_panel_hits.tsv", sep="\t"); pan["gene"] = pan.accession.map(PANEL)
assert pan.gene.notna().all()
pan = pan.merge(gff, on=["isolate","protein"], how="left")
pan.to_csv("out/panel_pfam_ncbifam_domain_hits.tsv", sep="\t", index=False)

# ---------- (2) mucin panel ----------
KEEP = ["GH2","GH16","GH20","GH27","GH29","GH33","GH35","GH36","GH42","GH84","GH85","GH89","GH95","GH98","GH101","GH109","GH110","GH112",
        "GH123","GH129","GH136","GH141","GH151","GH156","CBM32","CBM40","CBM47","CBM51","CBM67","CE20"]
GROUP = {"GH33":"Sialidase GH33","GH29":"Fucosidase GH29/GH95","GH95":"Fucosidase GH29/GH95","GH20":"Hexosaminidase GH20/GH84/GH85",
         "GH84":"Hexosaminidase GH20/GH84/GH85","GH85":"Hexosaminidase GH20/GH84/GH85","GH2":"Galactosidase GH2/GH35/GH42","GH35":"Galactosidase GH2/GH35/GH42",
         "GH42":"Galactosidase GH2/GH35/GH42","GH101":"Endo-a-GalNAc'ase GH101/GH129","GH129":"Endo-a-GalNAc'ase GH101/GH129","GH16":"GH16 (mostly non-mucin)",
         "GH109":"a-GalNAc'ase GH109 (Gfo/Idh fold, low specificity)"}
rawm = raw[raw.family.isin(KEEP)].drop_duplicates(["isolate","protein","family"])   # run_hmm.py logic: subset search, no cross-family overlap removal
resm = fam[fam.family.isin(KEEP)]
rows = []
for f in KEEP:
    for i in ISO:
        rows.append((i, "dbCAN", f, GROUP.get(f, "other mucin-subset family"), int(((rawm.isolate==i)&(rawm.family==f)).sum()), int(((resm.isolate==i)&(resm.family==f)).sum())))
PF_M = {"PF00884":"Sulfatase","PF03781":"FGE-sulfatase","PF13402":"Peptidase_M60","PF10462":"Peptidase_M66 (StcE)","PF02128":"Peptidase_M36","PF00026":"Asp (A1)","PF03629":"SASA","PF04131":"NanE"}
for acc, n in PF_M.items():
    for i in ISO:
        k = pan[(pan.isolate==i)&(pan.accession==acc)].protein.nunique()
        rows.append((i, "Pfam", f"{n} {acc}", {"PF00884":"Sulfatase + FGE","PF03781":"Sulfatase + FGE","PF13402":"Mucinase M60","PF10462":"Mucinase M66 StcE"}.get(acc, "other Pfam"), k, k))
muc = pd.DataFrame(rows, columns=["isolate","db","feature","group","n_proteins","n_proteins_after_cazome_overlap_resolution"])
muc.to_csv("out/mucin_panel_long.tsv", sep="\t", index=False)
mw = muc.pivot_table(index=["group","feature"], columns="isolate", values="n_proteins", aggfunc="sum").reindex(columns=ISO)
mw.to_csv("out/mucin_panel_matrix.tsv", sep="\t")
gw = muc.drop_duplicates(["isolate","feature"]).groupby(["group","isolate"]).n_proteins.sum().unstack().reindex(columns=ISO)
gw.to_csv("out/mucin_panel_group_matrix.tsv", sep="\t")
mp = rawm.merge(gff, on=["isolate","protein"], how="left")[["isolate","protein","contig","start","end","strand","family","model","i_evalue","score","hmm_cov","prokka_product"]]
mp.to_csv("out/mucin_panel_protein_hits.tsv", sep="\t", index=False)
gh16 = raw[raw.family=="GH16"].sort_values("score", ascending=False).drop_duplicates(["isolate","protein"])
gh16 = gh16.merge(gff, on=["isolate","protein"], how="left")[["isolate","protein","model","score","i_evalue","prokka_product"]].rename(columns={"model":"best_GH16_subfamily_model"})
print("GH16 subfamilies:\n", gh16.groupby(["isolate","best_GH16_subfamily_model"]).size())

# ---------- operon neighbourhoods ----------
CORE = {"sialic_acid": ["nanA","nanK","nanE"], "L-fucose": ["fucI","fucK","fucA"], "GlcNAc_GalNAc": ["nagA","nagB"]}
ACC = {"sialic_acid": ["TRAP_DctP","TRAP_DctM","TRAP_DctQ","MFS_1","Sugar_tr","ABC_SBP_1","ABC_SBP_8","ABC_permease_BPD"],
       "L-fucose": ["fucO","fucU","fucP","RbsD_FucU(generic)","fucI_domain"], "GlcNAc_GalNAc": ["agaS","aga_PTS","kbaY_gatY","kbaZ_gatZ","lacC","lacD","Glucosamine_iso(generic)"]}
gp = pan.copy(); gp.loc[gp.gene=="nagB-II","gene"] = "nagB"
specific = set(sum(CORE.values(), [])) | {"fucO","fucU","fucP","agaS","aga_PTS","kbaY_gatY","kbaZ_gatZ","lacC","lacD"}
pg = gp.sort_values("seq_score", ascending=False).groupby(["isolate","protein"]).agg(genes=("gene", lambda s: ",".join(sorted(set(s)))), best_model=("accession","first")).reset_index()
pg = pg.merge(gff, on=["isolate","protein"])
cazloc = fam.groupby(["isolate","protein"]).family.agg(lambda s: ",".join(sorted(set(s)))).rename("cazy").reset_index().merge(gff, on=["isolate","protein"])
op_rows, nb_rows = [], []
for path in CORE:
    want = set(CORE[path]) | set(ACC[path])
    for i in ISO:
        cand = pg[(pg.isolate==i) & pg.genes.map(lambda s: bool(set(s.split(",")) & want)).astype(bool)]
        anchors = cand[cand.genes.map(lambda s: bool(set(s.split(",")) & set(CORE[path]))).astype(bool)]
        present = {gname: int(anchors.genes.str.split(",").map(lambda s: gname in s).sum()) for gname in CORE[path]}
        best = []; 
        for a in anchors.itertuples():  # neighbourhood around every core-gene anchor, +-10 kb
            w = cand[(cand.contig==a.contig) & (cand.end >= a.start - 10000) & (cand.start <= a.end + 10000)]
            core_here = sorted(set(sum([x.split(",") for x in w.genes], [])) & set(CORE[path]))
            acc_here = sorted(set(sum([x.split(",") for x in w.genes], [])) & set(ACC[path]))
            cz = cazloc[(cazloc.isolate==i)&(cazloc.contig==a.contig)&(cazloc.end >= a.start-10000)&(cazloc.start <= a.end+10000)]
            nb_rows.append((path, i, a.protein, a.genes, a.contig, a.start, a.end, ",".join(core_here), ",".join(acc_here),
                            ";".join(f"{x.protein}:{x.genes}" for x in w.itertuples() if x.protein != a.protein),
                            ";".join(f"{x.protein}:{x.cazy}" for x in cz.itertuples())))
            best.append(len(core_here))
        nc = max(best) if best else 0
        call = ("complete core cluster" if nc == len(CORE[path]) else "partial cluster" if nc >= 2 else
                "isolated core gene(s) only" if nc == 1 else "absent")
        op_rows.append((i, path, *[present[x] for x in CORE[path]] + [None]*(3-len(CORE[path])), nc, len(CORE[path]), call))
op = pd.DataFrame(op_rows, columns=["isolate","pathway","core1_n","core2_n","core3_n","max_core_genes_colocated_10kb","n_core_genes","call"])
op["core_genes"] = op.pathway.map(lambda p: ",".join(CORE[p]))
op.to_csv("out/catabolic_operon_calls.tsv", sep="\t", index=False)
pd.DataFrame(nb_rows, columns=["pathway","isolate","anchor_protein","anchor_genes","contig","start","end","core_genes_within_10kb","accessory_within_10kb","neighbour_panel_genes","CAZymes_within_10kb"]).to_csv("out/catabolic_operon_neighbourhoods.tsv", sep="\t", index=False)
print(op.to_string())
pg[pg.genes.map(lambda s: bool(set(s.split(",")) & specific))].to_csv("out/operon_gene_hits_specific_models.tsv", sep="\t", index=False)

# sub-threshold diagnostic for core operon models (best hit at E<=1e-5 regardless of TC/GA)
abc = Alphabet.amino(); diag = []
models = []
for acc in ["TIGR00683","NF003164","NF003461","NF047821","NF002231","PF04131","TIGR01089","NF008220","TIGR02628","TIGR01086","TIGR02638","NF011949","TIGR00221","NF008371","TIGR00502","NF046059"]:
    with HMMFile(f"db/hmm/{acc}.hmm") as f: models += list(f)
for i in ISO:
    with SequenceFile(f"fg_inputs/{i}/{i}.faa", digital=True, alphabet=abc) as sf: seqs = sf.read_block()
    for top in pyhmmer.hmmsearch(models, seqs, cpus=2, E=1e-5):
        q = top.query; acc = q.accession.decode() if isinstance(q.accession, bytes) else q.accession
        cut = q.cutoffs.trusted if not acc.startswith("PF") else q.cutoffs.gathering
        h = top[0] if len(top) else None
        diag.append((i, acc.split(".")[0], PANEL[acc.split(".")[0]], cut[0] if cut else None,
                     (h.name.decode() if isinstance(h.name, bytes) else h.name) if h else "", round(h.score, 1) if h else None, h.evalue if h else None))
dg = pd.DataFrame(diag, columns=["isolate","accession","gene","cutoff_bits","best_hit_protein","best_hit_bits","best_hit_evalue"])
dg = dg.merge(gff[["isolate","protein","prokka_product"]], left_on=["isolate","best_hit_protein"], right_on=["isolate","protein"], how="left").drop(columns="protein")
dg["passes_cutoff"] = dg.best_hit_bits >= dg.cutoff_bits
dg.to_csv("out/operon_model_best_hits_diagnostic.tsv", sep="\t", index=False)

# ---------- (3) fungal cell-wall panel ----------
FW = {"GH18":"chitin","GH19":"chitin","GH23":"chitin? (mostly peptidoglycan lytic transglycosylase)","AA10":"chitin/cellulose LPMO","CBM5":"chitin-binding",
      "CBM12":"chitin-binding","CBM73":"chitin-binding","GH20":"chitobiase/hexosaminidase (noted)","GH46":"chitosan (extra)","GH75":"chitosan (extra)",
      "GH16":"beta-1,3/1,4-glucan","GH55":"beta-1,3-glucan","GH64":"beta-1,3-glucan","GH81":"beta-1,3-glucan","GH17":"beta-1,3-glucan",
      "GH128":"beta-1,3-glucan","GH132":"beta-1,3-glucan"}
ds = pd.concat([pd.read_csv(f"deepsig_out/{i}.deepsig.tsv", sep="\t").assign(isolate=i) for i in ISO])
tat = pan[pan.gene=="Tat_signal"]; tat = tat[tat.env_from <= 60].groupby(["isolate","protein"]).size().rename("tat_hmm").reset_index()
fwp = fam[fam.family.isin(FW)].copy()
arch = fam.groupby(["isolate","protein"]).family.agg(lambda s: "+".join(sorted(set(s)))).rename("all_cazy_families_on_protein").reset_index()
sub16 = gh16.set_index(["isolate","protein"]).best_GH16_subfamily_model
fwp["target"] = fwp.family.map(FW)
fwp["model_subfamily"] = [sub16.get((r.isolate, r.protein), r.model) if r.family=="GH16" else r.model for r in fwp.itertuples()]
fwp = fwp[["isolate","protein","family","model_subfamily","target","i_evalue","score","prot_len","contig","start","end","strand","prokka_product"]]
fwp = fwp.merge(arch, on=["isolate","protein"]).merge(ds, on=["isolate","protein"], how="left").merge(tat, on=["isolate","protein"], how="left")
fwp["tat_hmm"] = fwp.tat_hmm.fillna(0).astype(int) > 0
fwp["predicted_secreted_SP"] = fwp.deepsig_class == "SignalPeptide"
fwp.to_csv("out/fungal_cellwall_protein_hits.tsv", sep="\t", index=False)
fl = fwp.groupby(["isolate","family"]).agg(n_proteins=("protein","nunique"), n_with_signal_peptide=("predicted_secreted_SP","sum")).reset_index()
fl["target"] = fl.family.map(FW)
full = pd.MultiIndex.from_product([ISO, list(FW)], names=["isolate","family"]).to_frame(index=False)
fl = full.merge(fl, how="left").fillna({"n_proteins":0,"n_with_signal_peptide":0}); fl["target"] = fl.family.map(FW)
fl[["n_proteins","n_with_signal_peptide"]] = fl[["n_proteins","n_with_signal_peptide"]].astype(int)
fl.to_csv("out/fungal_cellwall_panel_long.tsv", sep="\t", index=False)
a = fl.pivot(index="family", columns="isolate", values="n_proteins").reindex(index=list(FW), columns=ISO)
b = fl.pivot(index="family", columns="isolate", values="n_with_signal_peptide").reindex(index=list(FW), columns=ISO)
pd.concat({"n_proteins": a, "n_with_signal_peptide": b}, axis=1).to_csv("out/fungal_cellwall_panel_matrix.tsv", sep="\t")
print(pd.concat({"n": a, "SP": b}, axis=1).to_string())
spsum = ds.groupby("isolate").deepsig_class.value_counts().unstack(); print(spsum)
spsum.to_csv("out/deepsig_class_counts_by_isolate.tsv", sep="\t")
ds.to_csv("out/deepsig_all_proteins.tsv", sep="\t", index=False)

# ---------- addenda ----------
# near-threshold operon candidates: best hit >= 90% of the model's TC/GA cutoff but below it (reported separately, provisional)
dg["frac_of_cutoff"] = (dg.best_hit_bits / dg.cutoff_bits).round(3)
near = dg[(~dg.passes_cutoff) & (dg.frac_of_cutoff >= 0.9)]
PATH_OF = {g: p for p, gl in CORE.items() for g in gl}; PATH_OF["nagB-II"] = "GlcNAc_GalNAc"
near = near.assign(pathway=near.gene.map(PATH_OF))
op["near_threshold_candidates(>=0.9xTC)"] = [";".join(f"{r.gene}:{r.best_hit_protein}({r.best_hit_bits}/{r.cutoff_bits:.0f} bits)"
                                              for r in near[(near.isolate==o.isolate)&(near.pathway==o.pathway)].itertuples()) for o in op.itertuples()]
op["note"] = ""
op.loc[(op.isolate=="BF3")&(op.pathway=="sialic_acid"), "note"] = ("nanE IHHBNBAB_02553 adjacent to nagB IHHBNBAB_02554 on NODE_36 with a ROK-family kinase "
    "(02548, Prokka 'Fructokinase', below nanK TC) and ABC sugar-transporter genes (02549-02551, 02556-02557); no nanA detected; partial/provisional")
op.to_csv("out/catabolic_operon_calls.tsv", sep="\t", index=False)
dg.to_csv("out/operon_model_best_hits_diagnostic.tsv", sep="\t", index=False)
# mucin per-protein table with Pfam rows + DeepSig/Tat
pm = pan[pan.accession.isin(PF_M)].drop_duplicates(["isolate","protein","accession"])[["isolate","protein","contig","start","end","strand","accession","model","i_evalue","dom_score","hmm_cov","prokka_product"]]
pm = pm.rename(columns={"accession":"family","dom_score":"score"}); pm["family"] = pm.family.map(lambda a: f"{PF_M[a]} {a}")
mp2 = pd.concat([mp, pm], ignore_index=True)
mp2["group"] = mp2.family.map(lambda f: GROUP.get(f, {"Sulfatase PF00884":"Sulfatase + FGE","FGE-sulfatase PF03781":"Sulfatase + FGE","Peptidase_M60 PF13402":"Mucinase M60","Peptidase_M66 (StcE) PF10462":"Mucinase M66 StcE"}.get(f, "other")))
mp2 = mp2.merge(ds, on=["isolate","protein"], how="left").merge(tat, on=["isolate","protein"], how="left")
mp2["tat_hmm"] = mp2.tat_hmm.fillna(0).astype(int) > 0
mp2.to_csv("out/mucin_panel_protein_hits.tsv", sep="\t", index=False)
print(near[["isolate","gene","best_hit_protein","frac_of_cutoff"]].to_string())
# flag Gfo/Idh/MocA-fold families (hits are mostly NAD-dependent oxidoreductases, e.g. inositol dehydrogenases)
GFO = {"GH109","GH179","GH188"}
var = var.reset_index(); var["gfo_idh_fold_low_specificity"] = var.family.isin(GFO)
var.to_csv("out/cazy_family_variability_ranked.tsv", sep="\t", index=False)
t1 = var.head(10).assign(list="top10_all_families")
t2 = var[~var.gfo_idh_fold_low_specificity].head(10).assign(list="top10_excluding_Gfo-Idh_fold")
pd.concat([t1, t2]).to_csv("out/cazy_top10_most_variable_families.tsv", sep="\t", index=False)
print(pd.concat([t1, t2])[["list","class","family","log2_max1_min1","range"]+ISO].to_string())
