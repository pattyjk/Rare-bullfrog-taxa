#!/usr/bin/env python3
"""Combine AMRFinderPlus calls with curated KEGG metal-resistance KOs; add genomic context; draw figure.

Usage: python resistance_tables.py AMR_DIR FG_INPUTS_DIR kofam_all7_detail.tsv.gz TABLES_DIR OUTDIR
  AMR_DIR      output of run_amrfinder.sh
  FG_INPUTS    BFn/BFn.{faa,gff} Prokka files
  TABLES_DIR   outputs/08_functional_genomics/tables (islands_BF1/7, isolate_specific_clusters_BF1/7)
Outputs: resistance_genes_BF1-BF7.xlsx, resistance_genes_BF1-BF7.png
"""
import sys, os, re, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib as mpl, matplotlib.pyplot as plt
amr_dir, fg, kofam, tabs, out = sys.argv[1:6]; os.makedirs(out, exist_ok=True)
ISO = [f'BF{i}' for i in range(1, 8)]

A = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(glob.glob(f'{amr_dir}/BF*.amrfinder.tsv'))])
det = pd.read_csv(kofam, sep='\t', low_memory=False)
det = det[det['KO'].astype(str).str.match(r'K\d{5}')].copy()
for c in ['thrshld', 'score']: det[c] = pd.to_numeric(det[c], errors='coerce')
det['iso'] = det['gene name'].str.split('__').str[0]; det['locus'] = det['gene name'].str.split('__').str[1]
det['sig'] = det['#'].astype(str).str.strip().eq('*')

cur = {  # curated metal-resistance KOs: efflux/detoxification/regulators only (no uptake, no B12, no proteases)
 'K00537': ('Arsenic', 'arsC (glutaredoxin)'), 'K03741': ('Arsenic', 'arsC (thioredoxin)'), 'K18701': ('Arsenic', 'arsC (mycothiol)'),
 'K03325': ('Arsenic', 'acr3 arsenite efflux'), 'K03893': ('Arsenic', 'arsB arsenite pump'), 'K01551': ('Arsenic', 'arsA ATPase*'),
 'K03892': ('Arsenic', 'arsR regulator'), 'K07755': ('Arsenic', 'arsM methyltransferase'), 'K11811': ('Arsenic', 'arsH'),
 'K08363': ('Mercury', 'merT'), 'K08365': ('Mercury', 'merR'),
 'K17686': ('Copper/silver', 'copA Cu+ ATPase'), 'K01533': ('Copper/silver', 'copB Cu2+ ATPase'), 'K07156': ('Copper/silver', 'copC'),
 'K07233': ('Copper/silver', 'copB (outer membrane)'), 'K07245': ('Copper/silver', 'copD'), 'K07787': ('Copper/silver', 'cusA'),
 'K07798': ('Copper/silver', 'cusB'), 'K07796': ('Copper/silver', 'cusC'), 'K07810': ('Copper/silver', 'cusF'), 'K06201': ('Copper/silver', 'cutC'),
 'K07213': ('Copper/silver', 'copZ chaperone'), 'K09796': ('Copper/silver', 'pcuC chaperone'), 'K19591': ('Copper/silver', 'cueR regulator'),
 'K21600': ('Copper/silver', 'csoR regulator'), 'K22041': ('Copper/silver', 'TetR Cu repressor'),
 'K01534': ('Zinc/cadmium/cobalt/nickel', 'zntA/cadA ATPase'), 'K15725': ('Zinc/cadmium/cobalt/nickel', 'czcC'),
 'K15726': ('Zinc/cadmium/cobalt/nickel', 'czcA'), 'K15727': ('Zinc/cadmium/cobalt/nickel', 'czcB'),
 'K16264': ('Zinc/cadmium/cobalt/nickel', 'czcD CDF efflux'), 'K24078': ('Zinc/cadmium/cobalt/nickel', 'cnrR'),
 'K21885': ('Zinc/cadmium/cobalt/nickel', 'cadC regulator'), 'K21886': ('Zinc/cadmium/cobalt/nickel', 'nmtR regulator'),
 'K22298': ('Zinc/cadmium/cobalt/nickel', 'smtB regulator'),
 'K07240': ('Chromate', 'chrA efflux'), 'K19784': ('Chromate', 'chrR reductase'),
 'K05791': ('Tellurite', 'terZ'), 'K05792': ('Tellurite', 'terA'), 'K05793': ('Tellurite', 'terB'), 'K05795': ('Tellurite', 'terD'),
 'K16868': ('Tellurite', 'tehB')}
k = det[det.sig & det.KO.isin(cur)].copy()
k['ratio'] = k.score / k.thrshld
k = k.sort_values('ratio', ascending=False).drop_duplicates(['iso', 'locus'])   # one call per protein (best KO)
k['category'] = k.KO.map(lambda x: cur[x][0]); k['gene'] = k.KO.map(lambda x: cur[x][1]); k['source'] = 'KEGG KofamScan'
a = A.rename(columns={'Name': 'iso', 'Protein id': 'locus', 'Element symbol': 'gene'}).copy()
a['category'] = np.where(a.Type == 'AMR', 'Antibiotic: ' + a.Class.str.title(), 'Arsenic'); a['source'] = 'AMRFinderPlus (HMM)'; a['KO'] = ''
bl = pd.read_csv(f'{amr_dir}/IHHBNBAB_01063_blastp.tsv', sep='\t', header=None, names=['q', 's', 'pid', 'len', 'qlen', 'slen', 'ev', 'bits']).iloc[0]
manual = pd.DataFrame([dict(iso='BF3', locus='IHHBNBAB_01063', gene=bl.s.split('|')[4], category='Antibiotic: Rifamycin', KO='',
                            source=f'manual BLASTP vs AMRProt ({bl.pid:.1f}% id, {bl.len}/{bl.slen} aa); AMRFinderPlus crashed on this protein')])
T = pd.concat([k[['iso', 'locus', 'KO', 'gene', 'category', 'source']],
               a[['iso', 'locus', 'KO', 'gene', 'category', 'source', '% Identity to reference', '% Coverage of reference']], manual], ignore_index=True)
T = T.drop_duplicates(['iso', 'locus'], keep='last')
T['class'] = np.where(T.category.str.startswith('Antibiotic'), 'Antibiotic resistance', 'Metal resistance')

rows = []
for i in ISO:
    for l in open(f'{fg}/{i}/{i}.gff'):
        if l.startswith('##FASTA'): break
        if l.startswith('#'): continue
        p = l.rstrip('\n').split('\t')
        if len(p) < 9 or p[2] != 'CDS': continue
        m = re.search(r'ID=([^;]+)', p[8]); pr = re.search(r'product=([^;]+)', p[8])
        rows.append((i, m.group(1), p[0], int(p[3]), int(p[4]), pr.group(1) if pr else ''))
GFF = pd.DataFrame(rows, columns=['iso', 'locus', 'contig', 'start', 'end', 'prokka_product'])
T = T.merge(GFF, on=['iso', 'locus'], how='left')
cov = lambda c: float(re.search(r'cov_([\d.]+)', c).group(1))
rel = {}
for i in ISO:
    cs = GFF[GFF.iso == i].contig.unique(); L = np.array([int(re.search(r'length_(\d+)', c).group(1)) for c in cs])
    cv = np.array([cov(c) for c in cs]); o = np.argsort(cv); med = cv[o][np.searchsorted(np.cumsum(L[o]), L.sum() / 2)]
    rel.update({(i, c): round(cov(c) / med, 2) for c in cs})
T['contig_cov_rel_median'] = [rel.get((i, c)) for i, c in zip(T.iso, T.contig)]
T['plasmid_like_contig(>=1.5x)'] = T.contig_cov_rel_median >= 1.5
I = {b: pd.read_csv(f'{tabs}/islands_{b}.tsv', sep='\t') for b in ['BF1', 'BF7']}
spec = set().union(*[set(pd.read_csv(f'{tabs}/isolate_specific_clusters_{b}.tsv', sep='\t').prokka_locus_tag.dropna()) for b in ['BF1', 'BF7']])
T['isolate_specific_gene'] = np.where(T.iso.isin(['BF1', 'BF7']), T.locus.isin(spec).map({True: 'yes', False: 'no'}), 'not assessed')
T['island'] = [';'.join(I[r.iso][(I[r.iso].contig == r.contig) & (I[r.iso].start <= r.end) & (I[r.iso].end >= r.start)].island_id)
               if r.iso in I else '' for r in T.itertuples()]

order = ['Antibiotic: Beta-Lactam', 'Antibiotic: Aminoglycoside', 'Antibiotic: Rifamycin', 'Antibiotic: Phenicol', 'Antibiotic: Macrolide',
         'Arsenic', 'Copper/silver', 'Zinc/cadmium/cobalt/nickel', 'Mercury', 'Chromate', 'Tellurite']
P = T.groupby(['category', 'iso']).size().unstack(fill_value=0).reindex(columns=ISO, fill_value=0).reindex(order).fillna(0).astype(int)
P['total'] = P.sum(axis=1)
nprot = pd.Series({i: sum(1 for l in open(f'{fg}/{i}/{i}.faa') if l.startswith('>')) for i in ISO})
print(P.to_string()); print('per isolate', T.groupby('iso').size().to_dict())
for b in ['BF1', 'BF7']:
    d = T[(T.iso == b) & (T['class'] == 'Metal resistance')]
    print(b, 'metal', len(d), 'in island', int((d.island != '').sum()), 'isolate-specific', int((d.isolate_specific_gene == 'yes').sum()))

outT = T[['iso', 'class', 'category', 'gene', 'locus', 'KO', 'source', '% Identity to reference', '% Coverage of reference', 'contig', 'start', 'end',
          'prokka_product', 'contig_cov_rel_median', 'plasmid_like_contig(>=1.5x)', 'island', 'isolate_specific_gene']].sort_values(['iso', 'class', 'category', 'gene'])
notes = ['Antibiotic resistance: NCBI AMRFinderPlus (see amrfinder_versions.txt; 4.2.7, db 2026-10-29.1 used), --plus, protein+nucleotide+Prokka GFF, no organism option. All calls are HMM family-level matches (45-77% identity); none is a close BLAST match to an acquired gene, so they are most likely intrinsic.',
         'BF3 protein IHHBNBAB_01063 crashes AMRFinderPlus 4.2.7 (amr_report HMM-domain parser); BF3 was run without it and the protein was called from BLASTP against AMRProt.fa.',
         'General multidrug efflux pumps are not reported by AMRFinderPlus and are not counted.',
         'Metal resistance: AMRFinderPlus STRESS/METAL calls plus curated KEGG KOs (KofamScan 1.3.0, above KO threshold); uptake transporters, cobalamin synthesis and metalloproteases excluded.',
         'Proteins hitting two curated KOs are kept once under the higher score/threshold KO.',
         'K01551 (arsA) also covers the GET3 tail-anchored protein ATPase and may not be arsenic-specific.',
         'contig_cov_rel_median: contig coverage / length-weighted median for that genome; >=1.5x = plasmid-like contig.',
         'island / isolate_specific_gene: from the BF1/BF7 comparative analysis (outputs/08); other isolates not assessed.',
         'Gene presence = capacity only; no MICs or metal tolerance were measured.']
with pd.ExcelWriter(f'{out}/resistance_genes_BF1-BF7.xlsx') as w:
    pd.DataFrame({'note': notes}).to_excel(w, sheet_name='methods', index=False)
    P.to_excel(w, sheet_name='counts_by_category')
    outT.to_excel(w, sheet_name='all_genes', index=False)
    A.to_excel(w, sheet_name='AMRFinderPlus_raw', index=False)

plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'axes.titlesize': 8})
gen = {'BF1': 'Sphingomonas', 'BF2': 'Xanthobact.', 'BF3': 'Kineococcus', 'BF4': 'Cryobacterium', 'BF5': 'Arthrobacter', 'BF6': 'Psychromicrob.', 'BF7': 'Rhodococcoides'}
fig, axs = plt.subplots(2, 1, figsize=(6.6, 5.6), gridspec_kw=dict(height_ratios=[1, 3.2], hspace=0.5))
ax = axs[0]
met = T[T['class'] == 'Metal resistance'].groupby('iso').size().reindex(ISO, fill_value=0) / nprot * 1000
ab = T[T['class'] == 'Antibiotic resistance'].groupby('iso').size().reindex(ISO, fill_value=0) / nprot * 1000
x = np.arange(7)
ax.bar(x, met.values, color='#117777', width=.6, label='Metal resistance')
ax.bar(x, ab.values, bottom=met.values, color='#771155', width=.6, label='Antibiotic resistance')
for xi, i in enumerate(ISO): ax.text(xi, met[i] + ab[i] + .5, f'n={int((T.iso == i).sum())}', ha='center', fontsize=6)
ax.set_xticks(x); ax.set_xticklabels(ISO); ax.set_ylabel('Genes per\n1,000 proteins'); ax.set_ylim(0, max(met + ab) * 1.3)
ax.legend(frameon=False, fontsize=6, loc='upper right', ncol=2)
ax.set_title('BF1 carries the largest resistance repertoire, mostly metal genes', loc='left')
ax.text(-0.12, 1.1, 'a', transform=ax.transAxes, fontsize=10, fontweight='bold')
ax = axs[1]; H = P[ISO]; v = H.values; vmax = v.max()
ax.imshow(np.log1p(v), cmap='Blues', aspect='auto', vmin=0, vmax=np.log1p(vmax) * 1.05)
for r in range(v.shape[0]):
    for c in range(v.shape[1]):
        ax.text(c, r, str(v[r, c]) if v[r, c] else '·', ha='center', va='center', fontsize=6,
                color='white' if np.log1p(v[r, c]) > np.log1p(vmax) * 0.6 else 'black')
ax.axhline(4.5, color='k', lw=0.8)
ax.set_yticks(range(len(order))); ax.set_yticklabels([o.replace('Antibiotic: ', '').replace('Beta-Lactam', 'β-lactam') for o in order])
ax.set_xticks(range(7))
ax.set_xticklabels([f'{i}\n{gen[i]}' if i == 'BF2' else f'{i}\n$\\it{{{gen[i].rstrip(".")}}}$' + ('.' if gen[i].endswith('.') else '') for i in ISO])
ax.text(-0.33, 2, 'Antibiotic', transform=ax.get_yaxis_transform(), rotation=90, va='center', ha='center', fontsize=7)
ax.text(-0.33, 8, 'Metal', transform=ax.get_yaxis_transform(), rotation=90, va='center', ha='center', fontsize=7)
ax.set_title('Resistance genes by class (proteins; antibiotic genes are family-level HMM matches)', loc='left')
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
ax.text(-0.12, 1.03, 'b', transform=ax.transAxes, fontsize=10, fontweight='bold')
fig.savefig(f'{out}/resistance_genes_BF1-BF7.png', dpi=300, bbox_inches='tight')
print('done')
