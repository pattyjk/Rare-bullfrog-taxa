#!/usr/bin/env python3
"""Dormancy-related genes in BF1-BF7 from KofamScan calls (plus Prokka names for universal stress proteins).

Usage: python dormancy_screen.py kofam_all7_detail.tsv.gz FG_INPUTS_DIR OUTDIR
  FG_INPUTS_DIR holds BFn/BFn.gff (Prokka)
Outputs: dormancy_genes_BF1-BF7.xlsx, dormancy_genes_BF1-BF7.png
"""
import sys, os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
kofam, fg, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
ISO = [f'BF{i}' for i in range(1, 8)]
det = pd.read_csv(kofam, sep='\t', low_memory=False)
det = det[det['KO'].astype(str).str.match(r'K\d{5}')].copy()
det['iso'] = det['gene name'].str.split('__').str[0]; det['sig'] = det['#'].astype(str).str.strip().eq('*')
D = [
 ('Endospore formation', 'spo0A master regulator', ['K07699']),
 ('Endospore formation', 'sporulation sigma (σE/F/G)', ['K03091']),
 ('Resuscitation from dormancy', 'Rpf (KEGG RpfA/B/E)', ['K21687', 'K21688', 'K21691']),
 ('DosR/DevSR dormancy regulon', 'DevS/DosS sensor', ['K07682']),
 ('DosR/DevSR dormancy regulon', 'DevR/DosR regulator', ['K07695']),
 ('DosR/DevSR dormancy regulon', 'small heat shock HSP20/α-crystallin (HspX-like)', ['K13993']),
 ('DosR/DevSR dormancy regulon', 'Tgs triacylglycerol synthase (WS/DGAT)', ['K00635']),
 ('Stringent response', 'Rel/SpoT (ppGpp synthase/hydrolase)', ['K01139']),
 ('Stringent response', 'small alarmone synthetase', ['K07816']),
 ('Ribosome hibernation', 'Hpf hibernation factor', ['K05808']),
 ('Ribosome hibernation', 'Rmf modulation factor', ['K03812']),
 ('Persistence (toxin–antitoxin)', 'toxins (MazF, RelE, ParE, HigB, HipA, SocB)', ['K07171', 'K06218', 'K19092', 'K07334', 'K07154', 'K28190']),
 ('Persistence (toxin–antitoxin)', 'antitoxins (VapB, RelB, ParD, HigA, CptB, SocA, TakA)',
  ['K18829', 'K28150', 'K28276', 'K28280', 'K18918', 'K28297', 'K07746', 'K21498', 'K09159', 'K28191', 'K28189']),
 ('Energy/carbon storage', 'PHB synthesis (phaC)', ['K03821']),
 ('Energy/carbon storage', 'glycogen (glgC)', ['K00975']),
 ('Energy/carbon storage', 'trehalose (otsA)', ['K00697']),
 ('Energy/carbon storage', 'trehalose (treZ, TreY/TreZ route)', ['K01236']),
 ('Energy/carbon storage', 'polyphosphate kinase (ppk1/ppk2)', ['K00937', 'K22468']),
 ('Stress survival', 'universal stress proteins (Usp)', None),   # Prokka product names
 ('Stress survival', 'Dps DNA protection', ['K04047']),
 ('Stress survival', 'cold-shock proteins', ['K03704']),
 ('Stress survival', 'diadenylate cyclase (c-di-AMP)', ['K07067'])]
prod = []
for i in ISO:
    for l in open(f'{fg}/{i}/{i}.gff'):
        if l.startswith('##FASTA'): break
        p = l.rstrip('\n').split('\t')
        if len(p) >= 9 and p[2] == 'CDS':
            m = re.search(r'product=([^;]+)', p[8]); prod.append((i, m.group(1) if m else ''))
prod = pd.DataFrame(prod, columns=['iso', 'product'])
rows = []
for s, g, ks in D:
    if ks is None:
        c = prod[prod['product'].str.contains(r'^Universal stress protein', case=False, regex=True)].groupby('iso').size()
        src, kk = 'Prokka product', ''
    else:
        c = det[det.sig & det.KO.isin(ks)].drop_duplicates(['iso', 'gene name']).groupby('iso').size(); src, kk = 'KEGG', ','.join(ks)
    r = {'strategy': s, 'gene': g, 'KOs': kk, 'source': src}; r.update(c.reindex(ISO, fill_value=0).to_dict()); rows.append(r)
DT = pd.DataFrame(rows)
print(DT[['gene'] + ISO].to_string(index=False))
notes = ['Counts = proteins above KEGG KO threshold (KofamScan 1.3.0) per genome, except universal stress proteins (Prokka product names).',
 'Endospore: spo0A (K07699) and sporulation sigma (K03091) are below threshold in all isolates (best 0.63-0.74). Scattered spore-associated KOs occur without a sporulation program and are not counted.',
 'Rpf: KEGG RpfA/B/E (K21687/K21688/K21691).',
 'DevS/DevR: at threshold in BF3-BF5; BF6/BF7 best scores 0.59-0.89. HSP20 family (K13993) includes HspX/Acr but also generic small heat-shock proteins.',
 'PHB: phaC (K03821) only BF1 and BF2, with phaA/phaB also present in both.',
 'Toxin-antitoxin: summed over the KOs listed; some may be addiction modules of mobile elements rather than dormancy systems.',
 'Gene presence = capacity only; dormancy was not measured.']
with pd.ExcelWriter(f'{out}/dormancy_genes_BF1-BF7.xlsx') as w:
    pd.DataFrame({'note': notes}).to_excel(w, sheet_name='methods', index=False)
    DT.to_excel(w, sheet_name='dormancy_genes', index=False)

plt.rcParams.update({'font.size': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 6})
gen = {'BF1': 'Sphingomonas', 'BF2': 'Xanthobact.', 'BF3': 'Kineococcus', 'BF4': 'Cryobacterium', 'BF5': 'Arthrobacter', 'BF6': 'Psychromicrob.', 'BF7': 'Rhodococcoides'}
v = DT[ISO].values
fig, ax = plt.subplots(figsize=(6.6, 6.2))
ax.imshow(np.log1p(v), cmap='Blues', aspect='auto', vmin=0, vmax=np.log1p(v.max()) * 1.05)
for r in range(v.shape[0]):
    for c in range(7):
        ax.text(c, r, str(v[r, c]) if v[r, c] else '·', ha='center', va='center', fontsize=6,
                color='white' if np.log1p(v[r, c]) > np.log1p(v.max()) * 0.6 else 'black')
ax.set_yticks(range(len(DT))); ax.set_yticklabels(DT.gene)
ax.set_xticks(range(7))
ax.set_xticklabels([f'{i}\n{gen[i]}' if i == 'BF2' else f'{i}\n$\\it{{{gen[i].rstrip(".")}}}$' + ('.' if gen[i].endswith('.') else '') for i in ISO])
for sp in ax.spines.values(): sp.set_visible(False)
ax.tick_params(length=0)
fig.canvas.draw(); rr = fig.canvas.get_renderer()
xa = ax.transAxes.inverted().transform((min(t.get_window_extent(rr).x0 for t in ax.get_yticklabels()) - 8, 0))[0]
for s, d in DT.groupby('strategy', sort=False):
    a, b = d.index.min(), d.index.max()
    if a > 0: ax.axhline(a - 0.5, color='k', lw=0.6)
    ax.text(xa, (a + b) / 2, s, transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=6, fontweight='bold')
ax.set_title('No isolate has an endospore program; Rpf and DosR genes occur only in the Actinobacteria (genes per genome)', loc='left', fontsize=8)
fig.savefig(f'{out}/dormancy_genes_BF1-BF7.png', dpi=300, bbox_inches='tight')
print('done')
