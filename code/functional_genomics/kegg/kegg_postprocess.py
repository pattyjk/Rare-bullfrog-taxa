#!/usr/bin/env python3
"""Post-process KofamScan detail-tsv (combined 7-proteome run, IDs prefixed '<ISOLATE>__').
Keeps only significant assignments ('*' = score >= KOfam adaptive threshold).
Outputs: per-protein KO table, KO long/wide count matrices, module completeness (via
module_completeness.py), curated summary tables. Run from kegg/ directory."""
import csv, re, os, sys, collections
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import module_completeness as mc
ISO = ['BF1','BF2','BF3','BF4','BF5','BF6','BF7']
TAXA = {'BF1':'Sphingomonas paucimobilis','BF2':'Xanthobacteraceae sp. nov.','BF3':'Kineococcus sp. nov.',
        'BF4':'Cryobacterium sp. nov.','BF5':'Arthrobacter sp. nov.','BF6':'Psychromicrobium sp. nov.','BF7':'Rhodococcoides fascians'}
os.makedirs('out', exist_ok=True)
# contig per protein from FAA headers
contig = {}
for i in ISO:
    for l in open(f'fg_inputs/{i}/{i}.faa'):
        if l.startswith('>'):
            p = l[1:].split(); m = re.search(r'contig=(\S+)', l)
            contig[f'{i}__{p[0]}'] = m.group(1) if m else ''
# 1. per-protein significant KO assignments
rows = []
for l in open('kofam_all7_detail.tsv'):
    if not l.startswith('*'): continue
    f = l.rstrip('\n').split('\t')
    _, gene, ko, thr, score, ev, desc = f[:7]
    iso, locus = gene.split('__', 1)
    rows.append(dict(isolate=iso, locus_tag=locus, contig=contig.get(gene, ''), KO=ko, threshold=float(thr),
                     score=float(score), evalue=float(ev), KO_definition=desc.strip('"')))
pp = pd.DataFrame(rows)
nko = pp.groupby(['isolate','locus_tag']).KO.transform('nunique'); pp['n_sig_KOs_for_protein'] = nko
pp = pp.sort_values(['isolate','locus_tag','score'], ascending=[True,True,False])
pp.to_csv('out/kofam_per_protein_KO.tsv', sep='\t', index=False)
# 2. KO x isolate counts (proteins per KO)
cnt = pp.groupby(['isolate','KO']).locus_tag.nunique().rename('n_proteins').reset_index()
defs = pp.drop_duplicates('KO').set_index('KO').KO_definition
cnt['KO_definition'] = cnt.KO.map(defs)
cnt.to_csv('out/KO_isolate_long.tsv', sep='\t', index=False)
wide = cnt.pivot(index='KO', columns='isolate', values='n_proteins').reindex(columns=ISO).fillna(0).astype(int)
wide.insert(0, 'KO_definition', wide.index.map(defs))
wide['n_isolates_present'] = (wide[ISO] > 0).sum(axis=1)
wide.to_csv('out/KO_isolate_count_matrix.tsv', sep='\t')
pres = cnt[['isolate','KO']]; pres.to_csv('out/_ko_presence.tsv', sep='\t', index=False)
# 3. modules
mc.main('modules_raw', 'out/_ko_presence.tsv', 'out/module_completeness')
ml = pd.read_csv('out/module_completeness_long.tsv', sep='\t')
# --- near-threshold diagnostics (NOT used for primary calls) ---
RATIO, EMAX = 0.5, 1e-10   # relaxed sensitivity: non-significant hit with score >= 0.5 x adaptive threshold and E <= 1e-10
dd = pd.read_csv('kofam_all7_detail.tsv', sep='\t', skiprows=[1], header=0, quoting=3)
dd.columns = ['sig','gene','KO','thr','score','E','def']
dd['thr'] = pd.to_numeric(dd.thr, errors='coerce'); dd = dd.dropna(subset=['thr'])
dd['score'] = pd.to_numeric(dd.score); dd['E'] = pd.to_numeric(dd.E)
dd['isolate'] = dd.gene.str.split('__').str[0]; dd['locus_tag'] = dd.gene.str.split('__').str[1]
dd['ratio'] = dd.score / dd.thr
ns = dd[dd.sig.astype(str).str.strip() != '*']
best_ns = ns.sort_values('score', ascending=False).drop_duplicates(['isolate','KO']).set_index(['isolate','KO'])
rel = ns[(ns.ratio >= RATIO) & (ns.E <= EMAX)][['isolate','KO']].drop_duplicates()
pd.concat([pres, rel]).drop_duplicates().to_csv('out/_ko_presence_relaxed.tsv', sep='\t', index=False)
mc.main('modules_raw', 'out/_ko_presence_relaxed.tsv', 'out/_relaxed')
rlx = pd.read_csv('out/_relaxed_matrix.tsv', sep='\t').set_index('module')
mw = pd.read_csv('out/module_completeness_matrix.tsv', sep='\t')
cls = mw['class'].str.split('; ', expand=True)
mw.insert(3, 'module_type', cls[0]); mw.insert(4, 'class_l2', cls[1]); mw.insert(5, 'class_l3', cls[2])
for i in ISO: pass
mw['n_isolates_complete'] = (mw[ISO] >= 0.9999).sum(axis=1); mw['n_isolates_ge075'] = (mw[ISO] >= 0.75).sum(axis=1)
mw.to_csv('out/module_completeness_matrix.tsv', sep='\t', index=False)
# 4. curated module summary
CUR = [
 ('Vitamin/cofactor','B1 thiamine',['M00127','M00895','M00898','M00899']),
 ('Vitamin/cofactor','B2 riboflavin',['M00125']),
 ('Vitamin/cofactor','B3 NAD de novo',['M00115','M00912']),
 ('Vitamin/cofactor','B5 pantothenate',['M00119','M00913']),
 ('Vitamin/cofactor','CoA from pantothenate',['M00120']),
 ('Vitamin/cofactor','B6 pyridoxal-P',['M00124','M00916']),
 ('Vitamin/cofactor','B7 biotin',['M00123','M00573','M00577','M00950']),
 ('Vitamin/cofactor','B9 tetrahydrofolate',['M00126','M00840','M00841']),
 ('Vitamin/cofactor','B12 cobalamin (corrin ring)',['M00924','M00925']),
 ('Vitamin/cofactor','B12 cobalamin (cobinamide => cobalamin)',['M00122']),
 ('Vitamin/cofactor','Heme',['M00121','M00926']),
 ('Vitamin/cofactor','Siroheme',['M00846']),
 ('Vitamin/cofactor','Menaquinone',['M00116','M00930','M00931']),
 ('Vitamin/cofactor','Ubiquinone',['M00117']),
 ('Vitamin/cofactor','Molybdenum cofactor',['M00880']),
 ('Amino acid','Serine',['M00020']),
 ('Amino acid','Cysteine',['M00021','M00338','M00609']),
 ('Amino acid','Threonine',['M00018']),
 ('Amino acid','Methionine',['M00017']),
 ('Amino acid','Lysine',['M00016','M00525','M00526','M00527','M00030','M00031']),
 ('Amino acid','Valine (+isoleucine from 2-oxobutanoate)',['M00019']),
 ('Amino acid','Isoleucine (2-oxobutanoate supply)',['M00535','M00570']),
 ('Amino acid','Leucine',['M00432']),
 ('Amino acid','Proline',['M00015']),
 ('Amino acid','Ornithine',['M00028','M00763']),
 ('Amino acid','Arginine',['M00844','M00845']),
 ('Amino acid','Histidine',['M00026']),
 ('Amino acid','Chorismate (shikimate)',['M00022']),
 ('Amino acid','Tryptophan',['M00023']),
 ('Amino acid','Phenylalanine',['M00024']),
 ('Amino acid','Tyrosine',['M00025']),
 ('Nitrogen','Assimilatory nitrate reduction',['M00531']),
 ('Nitrogen','Dissimilatory nitrate reduction (DNRA)',['M00530']),
 ('Nitrogen','Denitrification',['M00529']),
 ('Nitrogen','Nitrification',['M00528','M00804']),
 ('Nitrogen','Nitrogen fixation',['M00175']),
 ('Sulfur','Assimilatory sulfate reduction',['M00176']),
 ('Sulfur','Dissimilatory sulfate reduction',['M00596']),
 ('Sulfur','SOX thiosulfate oxidation',['M00595']),
 ('Central carbon','Glycolysis (EMP)',['M00001']),
 ('Central carbon','Glycolysis core (C3)',['M00002']),
 ('Central carbon','Gluconeogenesis',['M00003']),
 ('Central carbon','Pyruvate oxidation',['M00307']),
 ('Central carbon','TCA cycle',['M00009']),
 ('Central carbon','Entner-Doudoroff',['M00008']),
 ('Central carbon','Semi-phosphorylative ED',['M00308','M00633']),
 ('Central carbon','Glyoxylate cycle',['M00012']),
 ('Central carbon','Pentose phosphate (oxidative)',['M00006']),
 ('Central carbon','Pentose phosphate (non-oxidative)',['M00007']),
]
mwi = mw.set_index('module')
out = []
for cat, feat, mods in CUR:
    for i in ISO:
        vals = {m: mwi.loc[m, i] for m in mods if m in mwi.index}
        best = max(vals, key=vals.get)
        out.append(dict(category=cat, feature=feat, isolate=i, taxon=TAXA[i], best_module=best,
                        best_module_name=mwi.loc[best,'name'], completeness=round(vals[best],3),
                        complete=int(vals[best] >= 0.9999), ge075=int(vals[best] >= 0.75),
                        lt050_possible_auxotrophy=int(vals[best] < 0.5) if cat=='Amino acid' else '',
                        all_modules=';'.join(f'{m}={v:.2f}' for m, v in vals.items()),
                        sensitivity_completeness_relaxed=round(max(rlx.loc[m, i] for m in vals), 3)))
cs = pd.DataFrame(out)
cs.to_csv('out/curated_modules_long.tsv', sep='\t', index=False)
cmods = sorted({m for _, _, ms in CUR for m in ms})
mr = []
for r in ml[ml.module.isin(cmods)].itertuples():
    if not isinstance(r.missing_steps, str): continue
    for st in r.missing_steps.split(';'):
        for k in re.findall(r'K\d{5}', st):
            b = best_ns.loc[(r.isolate, k)] if (r.isolate, k) in best_ns.index else None
            mr.append(dict(isolate=r.isolate, module=r.module, module_name=r.name, missing_step=st, KO=k,
                           best_subthreshold_locus='' if b is None else b.locus_tag, score='' if b is None else b.score,
                           adaptive_threshold='' if b is None else b.thr, score_to_threshold='' if b is None else round(b.ratio, 3),
                           evalue='' if b is None else b.E))
pd.DataFrame(mr).to_csv('out/curated_missing_steps_subthreshold_hits.tsv', sep='\t', index=False)
print('relaxed-only KO assignments:', len(rel.merge(pres, how='left', indicator=True).query('_merge=="left_only"')))
# KO-level marker sets
kp = {(r.isolate, r.KO): r.n_proteins for r in cnt.itertuples()}
def n(i, k): return kp.get((i, k), 0)
KOSETS = [
 ('Nitrogen','Urease alpha ureC',['K01428']),('Nitrogen','Urease beta ureB',['K01429']),('Nitrogen','Urease gamma ureA',['K01430']),
 ('Nitrogen','Urease fused gamma-beta ureAB',['K14048']),
 ('Nitrogen','Urease accessory ureD/ureE/ureF/ureG',['K03190','K03187','K03188','K03189']),
 ('Nitrogen','Urea ABC transporter urtABCDE',['K11959','K11960','K11961','K11962','K11963']),
 ('Nitrogen','Urea carboxylase / allophanate hydrolase',['K01941','K01457']),
 ('Nitrogen','Glutamine synthetase glnA',['K01915']),
 ('Nitrogen','Glutamate synthase GOGAT gltB/gltD (NADPH)',['K00265','K00266']),
 ('Nitrogen','Glutamate synthase GOGAT (ferredoxin) / glutamate dehydrogenase',['K00284','K00261','K00262','K15371']),
 ('Nitrogen','Ammonium transporter amtB',['K03320']),
 ('Nitrogen','Assimilatory nitrate reductase narB/nasA/nasC',['K00367','K00372','K00360']),
 ('Nitrogen','Assimilatory nitrite reductase nirBD (NADH)',['K00362','K00363']),
 ('Nitrogen','Ferredoxin-nitrite reductase nirA',['K00366']),
 ('Nitrogen','Respiratory nitrate reductase narGHI',['K00370','K00371','K00374']),
 ('Nitrogen','Periplasmic nitrate reductase napAB',['K02567','K02568']),
 ('Nitrogen','NO-forming nitrite reductase nirK/nirS',['K00368','K15864']),
 ('Nitrogen','Nitric oxide reductase norBC',['K04561','K02305']),
 ('Nitrogen','Nitrous oxide reductase nosZ',['K00376']),
 ('Nitrogen','Nitrogenase nifHDK',['K02588','K02586','K02591']),
 ('Nitrogen','Nitrate/nitrite transporters narK/nrtABC',['K02575','K15576','K15577','K15578']),
 ('Nitrogen','Flavohemoprotein NO dioxygenase hmp',['K05916']),
]
ko_out = []
for cat, feat, kos in KOSETS:
    for i in ISO:
        c = [n(i, k) for k in kos]
        ko_out.append(dict(category=cat, feature=feat, isolate=i, taxon=TAXA[i], KOs=','.join(kos),
                           KOs_present=','.join(k for k, x in zip(kos, c) if x), n_KOs_present=sum(1 for x in c if x),
                           n_KOs_total=len(kos), n_proteins=sum(c)))
kdf = pd.DataFrame(ko_out)
for i in ISO:  # derived urease call
    s = {k: n(i, k) for k in ['K01428','K01429','K01430','K14048']}
    struct = s['K01428'] > 0 and ((s['K01429'] > 0 and s['K01430'] > 0) or s['K14048'] > 0)
    kdf.loc[(kdf.isolate==i)&(kdf.feature.str.startswith('Urease alpha')),'derived_call'] = 'urease structural genes complete' if struct else ('partial urease' if any(s.values()) else 'no urease structural genes')
kdf.to_csv('out/curated_KO_marker_sets_long.tsv', sep='\t', index=False)
# Pathway KO sets: steroid + aromatic degradation + TCS
lk = pd.read_csv('kegg_pathway_ko_links.tsv', sep='\t')
pk = []
for (pw, pn), g in lk.groupby(['pathway','pathway_name']):
    kos = sorted(set(g.KO))
    for i in ISO:
        pres_k = [k for k in kos if n(i, k)]
        pk.append(dict(pathway=pw, pathway_name=pn, isolate=i, taxon=TAXA[i], n_KOs_in_map=len(kos), n_KOs_present=len(pres_k),
                       n_proteins=sum(n(i,k) for k in pres_k), KOs_present=','.join(pres_k)))
pkd = pd.DataFrame(pk); pkd.to_csv('out/pathway_KO_counts_long.tsv', sep='\t', index=False)
st = lk[lk.pathway=='map00984'].KO.unique()
stw = wide.loc[wide.index.isin(st), ['KO_definition']+ISO+['n_isolates_present']]
stw.to_csv('out/steroid_degradation_map00984_KO_matrix.tsv', sep='\t')
TCS = [('PhoR/PhoB','K07636','K07657','phosphate limitation'),('PhoQ/PhoP','K07637','K07660','Mg2+ / cationic antimicrobial peptides / acidic pH'),
 ('BasS/BasR (PmrB/PmrA)','K07643','K07771','Fe3+ / acidic pH; LPS modification'),('QseC/QseB','K07645','K07666','host catecholamines / AI-3 (enteric literature)'),
 ('CpxA/CpxR','K07640','K07662','envelope stress'),('EnvZ/OmpR','K07638','K07659','osmolarity'),('KdpD/KdpE','K07646','K07667','K+ limitation / osmotic stress'),
 ('BaeS/BaeR','K07642','K07664','envelope stress / efflux'),('CusS/CusR','K07644','K07665','copper/silver'),('RstB/RstA','K07639','K07661','acid / PhoPQ network'),
 ('MtrB/MtrA','K07654','K07670','actinobacterial cell envelope / osmotic'),('SenX3/RegX3','K07768','K07776','phosphate (Actinobacteria)'),
 ('PrrB/PrrA','K07655','K07671','host-intracellular / N limitation (mycobacteria)'),('DevS/DevR (DosS/DosR)','K07682','K07695','hypoxia / NO / CO'),
 ('ChvG/ChvI','K14980','K14981','acidic pH; plant-host interaction (Alphaproteobacteria)'),('NarX/NarL','K07673','K07684','nitrate/nitrite'),
 ('NtrB/NtrC','K07708','K07712','nitrogen limitation'),('FixL/FixJ','K14986','K14987','oxygen'),('RegB/RegA','K15011','K15012','redox'),
 ('LiaS/LiaR','K11617','K11618','cell-envelope stress / antimicrobial peptides'),('VraS/VraR','K07681','K07694','cell-wall-active antibiotics'),
 ('VanS/VanR','K18345','K18344','glycopeptide'),('GacS/GacA (BarA/UvrY)','K07678','K07689','secondary-metabolite / biocontrol regulation'),
 ('CheA/CheY','K03407','K03413','chemotaxis'),('PilS/PilR','K02668','K02667','pilin expression / surface')]
tc = []
for nm_, hk, rr, sig in TCS:
    for i in ISO:
        tc.append(dict(system=nm_, sensor_KO=hk, regulator_KO=rr, signal_reported=sig, isolate=i, taxon=TAXA[i],
                       sensor_n=n(i,hk), regulator_n=n(i,rr), pair_present=int(n(i,hk)>0 and n(i,rr)>0)))
tcd = pd.DataFrame(tc); tcd.to_csv('out/two_component_host_sensing_long.tsv', sep='\t', index=False)
# combined curated summary (long + wide)
a = cs.assign(value=cs.completeness, value_type='module completeness (best alternative)', detail=cs.best_module+' '+cs.best_module_name)
b = kdf.assign(value=kdf.n_KOs_present/kdf.n_KOs_total, value_type='fraction of marker KOs present', detail=kdf.KOs_present)
c = pkd[pkd.pathway!='map02020'].assign(category='Aromatic/steroid degradation (KO count)', feature=lambda d: d.pathway+' '+d.pathway_name,
       value=lambda d: d.n_KOs_present, value_type='n KOs present in map', detail=lambda d: d.KOs_present)
d = tcd.assign(category='Two-component (host sensing)', feature=tcd.system+' ('+tcd.signal_reported+')', value=tcd.pair_present,
       value_type='HK+RR pair present (1/0)', detail=tcd.sensor_KO+':'+tcd.sensor_n.astype(str)+','+tcd.regulator_KO+':'+tcd.regulator_n.astype(str))
cols = ['category','feature','isolate','taxon','value','value_type','detail']
summ = pd.concat([a[cols+['complete','ge075','lt050_possible_auxotrophy','sensitivity_completeness_relaxed']], b[cols], c[cols], d[cols]], ignore_index=True)
summ.to_csv('out/curated_summary_long.tsv', sep='\t', index=False)
sw = summ.pivot_table(index=['category','feature','value_type'], columns='isolate', values='value', aggfunc='first', sort=False).reindex(columns=ISO)
sw.round(3).to_csv('out/curated_summary_wide.tsv', sep='\t')
print('proteins with sig KO per isolate:', pp.groupby('isolate').locus_tag.nunique().to_dict())
print('distinct KOs per isolate:', (wide[ISO]>0).sum().to_dict())
