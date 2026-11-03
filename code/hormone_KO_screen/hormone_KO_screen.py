#!/usr/bin/env python3
"""Steroid (Mohamed et al. 2026) and broader hormone-related KO screen of BF1-BF7.

Input : kofam_all7_detail.tsv.gz  (KofamScan 1.3.0 --format detail-tsv, all 7 genomes
        concatenated; gene name = '<isolate>__<locus_tag>'; '*' in column '#' = above KO threshold)
Output: steroid_genes_Mohamed2026_BF1-BF7.tsv, hormone_related_KOs_BF1-BF7.tsv
Usage : python hormone_KO_screen.py kofam_all7_detail.tsv.gz outdir/
"""
import sys, re, os
import numpy as np, pandas as pd

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
det = pd.read_csv(src, sep='\t', low_memory=False)
det = det[det['KO'].astype(str).str.match(r'K\d{5}')].copy()
for c in ['thrshld', 'score']:
    det[c] = pd.to_numeric(det[c], errors='coerce')
det['iso'] = det['gene name'].str.split('__').str[0]
det['sig'] = det['#'].astype(str).str.strip().eq('*')
ISO = [f'BF{i}' for i in range(1, 8)]
defs = det.groupby('KO')['KO definition'].first()

# 1. Mohamed et al. 2026 (Nat Microbiol, doi 10.1038/s41564-026-02494-1) steroid gene set
G = [('choD', 'K03333'), ('17βHSD', 'K05296'), ('tesI', 'K16051'), ('ksi', 'K01822'),
     ('kstD', 'K05898'), ('kshA', 'K15982'), ('kshB', 'K15983'), ('hsaA', 'K16047'),
     ('hsaC', 'K16049'), ('hsaD', 'K16050'), ('fadD3', 'K18687')]
rows = []
for g, k in G:
    d = det[det.KO == k]
    r = {'gene': g, 'KO': k, 'KO_definition': defs.get(k, '')}
    for i in ISO:
        di = d[d.iso == i]
        r[i] = int(di.sig.sum())
        r[i + '_best_ratio'] = round((di.score / di.thrshld).max(), 2) if len(di) else np.nan
    rows.append(r)
pd.DataFrame(rows).to_csv(f'{out}/steroid_genes_Mohamed2026_BF1-BF7.tsv', sep='\t', index=False)

# 2. Broader hormone/signalling screen: regex over KO definitions, significant hits only
pat = (r'indole-3-pyruvate decarboxylase|tryptophan 2-monooxygenase|indoleacetamide hydrolase|'
       r'aminocyclopropane-1-carboxylate deaminase|glutamate decarboxylase|aromatic-L-amino-acid decarboxylase|'
       r'tyrosine decarboxylase|histidine decarboxylase|cytokinin|tRNA dimethylallyltransferase|'
       r'adenylate dimethylallyltransferase|choloylglycine hydrolase|bile acid|desmolase|hydroxysteroid|'
       r'estradiol|testosterone|androst|ent-kaurene|copalyl|gibberellin|melatonin|serotonin|arylalkylamine|'
       r'tryptophan decarboxylase|phenylethylamine|dopamine|tyramine|cortisol|corticoster|7alpha-hydroxy|'
       r'estrone|indole-3-acetaldehyde|indoleacetate|IAA|tryptophan 5-monooxygenase|thyro|deiodinase|abscisic')
H = [k for k, d in defs.items() if re.search(pat, str(d), re.I)]
print('hormone-pattern KOs:', len(H))
sg = det[det.sig & det.KO.isin(H)]
cnt = sg.groupby(['KO', 'iso']).size().unstack(fill_value=0).reindex(columns=ISO, fill_value=0)
cnt.insert(0, 'def', [defs[k][:70] for k in cnt.index])
cnt = cnt.drop(index=['K24696'], errors='ignore')   # mycofactocin precursor decarboxylase: not a hormone gene
cnt.reset_index().to_csv(f'{out}/hormone_related_KOs_BF1-BF7.tsv', sep='\t', index=False)
print(cnt.to_string())
