# session cell 109 | 2026-10-17 21:33 | env: python | language: python
# files written: []


fix={"gephryomantis":"gephyromantis","ischconemia":"ischnocnema","hypsiboas":"boana"}
fam_manual={"salamandra":"salamandridae","gephyromantis":"mantellidae","boana":"hylidae","brachycephalus":"brachycephalidae","litoria":"hylidae","espadarana":"centrolenidae",
 "silverstoneia":"dendrobatidae","hyalinobatrachium":"centrolenidae","anodontohyla":"microhylidae","hemidactylium":"plethodontidae","pristimantis":"craugastoridae","stumffia":"microhylidae",
 "ischnocnema":"brachycephalidae","bufo":"bufonidae","haddadus":"craugastoridae","dendrobates":"dendrobatidae","leiopelma":"leiopelmatidae","plethodonthyla":"microhylidae","allobates":"aromobatidae",
 "hylodes":"hylodidae","paradoxohyla":"microhylidae","telmatobius":"telmatobiidae","psychrophrynella":"craugastoridae","adenomera":"leptodactylidae","bryophryne":"craugastoridae",
 "gastrotheca":"hemiphractidae","oreobates":"craugastoridae","osteocephalus":"hylidae"}
m["hg"]=m.hg.replace(fix); m["hsp"]=m.hg+" "+m.hsp.str.split().str[1].fillna("")
m["hf"]=m.hg.map(lambda g:g2f.get(g,fam_manual.get(g)))
sp_set=set(mpg.hsp); g_set=set(mpg.hg); f_set=set(mpg.hf)
def level(r):
    if r.hsp in sp_set and not r.hsp.endswith((" sp"," aff"," cf")): return "same host species"
    if r.hg in g_set: return "same host genus"
    if r.hf in f_set: return "same host family"
    return "no relative in meta-analysis"
m["home_level"]=m.apply(level,axis=1)
mg=m[(m.sOTUs!="")]
print(mg.home_level.value_counts().to_dict())
# per-sample rel abundance of all needed sOTUs
need=set(";".join(mg.sOTUs).split(";")); miss=need-set(RM0.index); print("need",len(need),"missing from RM0",len(miss))
fn="meta/otu/amphib_only_otu_table.txt"
with open(fn) as f: samples=f.readline().rstrip("\n").split("\t")
ns=len(samples); colidx=[samples.index(s) for s in gsamples]
tot=pd.read_pickle("meta/totals.pkl").seq_reads if False else None
