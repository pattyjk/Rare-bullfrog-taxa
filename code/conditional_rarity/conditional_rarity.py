#!/usr/bin/env python3
"""Consistently vs conditionally rare sOTUs in the amphibian skin meta-analysis.

Inputs (user's own data, not bundled):
  amphib_only_otu_table.zip   (amphib_only_otu_table.txt; 90-bp deblur V4 sOTUs as row IDs)
  amphib_only_map.txt         (sample metadata)
  meta_analysis_taxonomy.txt  (BLAST top hit per sOTU vs NCBI 16S RefSeq)
  AmphiBac_meta_abundance.xlsx, sheet your_isolates (genome-extracted V4 of the isolates)
Outputs: conditional_rarity_meta.xlsx, conditional_rarity_meta.png, sOTU_classes_all.tsv.gz
Usage: python conditional_rarity.py OTU.zip MAP.txt TAX.txt ISOLATES.xlsx OUTDIR
"""
import sys, os, re, zipfile, itertools, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib as mpl, matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
otu_zip, map_f, tax_f, iso_xlsx, out = sys.argv[1:6]
os.makedirs(out, exist_ok=True)

# ---- OTU table: keep well-formed rows (90-bp ACGT IDs), drop all-zero rows, merge duplicate IDs
z = zipfile.ZipFile(otu_zip)
lines = z.open('amphib_only_otu_table.txt').read().decode().split('\n')
cols = lines[0].rstrip('\n').split('\t'); n = len(cols)
ids, rows = [], []
for l in lines[1:]:
    if not l.strip(): continue
    p = l.rstrip('\r').split('\t')
    if len(p) != n + 1 or not re.fullmatch(r'[ACGTN]+', p[0]): continue
    ids.append(p[0]); rows.append(np.array(p[1:], dtype=float))
X = pd.DataFrame(np.vstack(rows), index=ids, columns=cols)
X = X[X.sum(axis=1) > 0].groupby(level=0).sum()
tot = X.sum(axis=0); X = X[tot[tot >= 1000].index]            # samples with >=1,000 reads
mp = pd.read_csv(map_f, sep='\t', low_memory=False).set_index('#SampleID').reindex(X.columns)
ms = mp[mp.microbial_habitat_type == 'external'].copy()       # skin samples only
R = X[ms.index] / X[ms.index].sum(axis=0)
print('sOTUs', X.shape[0], 'skin samples', R.shape[1])

# ---- units = host species x dataset (DOI|country) x life stage, >=15 samples
ms['ds'] = ms.doi.str.replace('DOI: ', '') + '|' + ms.country
ms['unit'] = ms.host_scientific_name + ' | ' + ms.ds + ' | ' + ms.lifestage
u = ms.groupby('unit').size(); u = u[u >= 15]; units = u.index.tolist()
info = pd.DataFrame([k.split(' | ') for k in units], index=units, columns=['host', 'ds', 'stage'])
M = pd.DataFrame({k: R[ms.index[ms.unit == k]].mean(axis=1) for k in units})
P = pd.DataFrame({k: (R[ms.index[ms.unit == k]] > 0).mean(axis=1) for k in units})
MX = pd.DataFrame({k: R[ms.index[ms.unit == k]].max(axis=1) for k in units})
keep = (M > 0).any(axis=1); M, P, MX = M[keep], P[keep], MX[keep]
C = pd.DataFrame(np.select([M >= 0.01, M >= 0.001, M > 0], ['A', 'I', 'R'], '0'), index=M.index, columns=units)
print('units', len(units), 'sOTUs in units', len(M))

# ---- pairwise: % of unit a's abundant sOTUs that are abundant / rare / absent in unit b
def differs(a, b):
    ia, ib = info.loc[a], info.loc[b]; d = []
    if ia.ds.split('|')[0] != ib.ds.split('|')[0]: d.append('study')
    if ia.ds.split('|')[1] != ib.ds.split('|')[1]: d.append('country')
    if ia.stage != ib.stage: d.append('stage')
    return '+'.join(d) or 'none'
rows = []
for a, b in itertools.permutations(units, 2):
    ab = C[a] == 'A'
    if ab.sum() == 0: continue
    cb = C.loc[ab, b]; ia, ib = info.loc[a], info.loc[b]
    if ia.host == ib.host: typ = 'same host, different context'
    elif ia.ds == ib.ds and ia.stage == ib.stage: typ = 'different host, same dataset/country/stage'
    else: typ = 'different host, different context'
    dif = differs(a, b) if ia.host == ib.host else ('host only' if typ.startswith('different host, same') else 'host+' + differs(a, b))
    rows.append(dict(a=a, b=b, type=typ, differs=dif, n_abund=int(ab.sum()),
                     pct_abund_in_b=(cb == 'A').mean() * 100, pct_rare_in_b=(cb == 'R').mean() * 100,
                     pct_absent_in_b=(cb == '0').mean() * 100))
PW = pd.DataFrame(rows); PW['rare_or_absent'] = PW.pct_rare_in_b + PW.pct_absent_in_b
print(PW.groupby('differs').rare_or_absent.agg(['size', 'median']).round(1).to_string())

# ---- per-sOTU classes
nU = len(units)
S = pd.DataFrame({'n_units_detected': (C != '0').sum(axis=1), 'n_abundant': (C == 'A').sum(axis=1),
                  'n_intermediate': (C == 'I').sum(axis=1), 'n_rare': (C == 'R').sum(axis=1)})
S['max_unit_mean_pct'] = M.max(axis=1) * 100
S['median_unit_prev_when_detected'] = P.where(P > 0).median(axis=1).round(3)
S['max_sample_pct'] = MX.max(axis=1) * 100
pairs_host = [(a, b) for a, b in itertools.permutations(units, 2) if info.loc[a, 'host'] != info.loc[b, 'host']
              and info.loc[a, 'ds'] == info.loc[b, 'ds'] and info.loc[a, 'stage'] == info.loc[b, 'stage']]
pairs_same = [(a, b) for a, b in itertools.permutations(units, 2) if info.loc[a, 'host'] == info.loc[b, 'host']]
def ev(pairs):
    o = pd.Series(False, index=C.index)
    for a, b in pairs: o |= (C[a] == 'A') & (C[b] == 'R')
    return o
S['host_contrast_AvsR'] = ev(pairs_host); S['same_host_contrast_AvsR'] = ev(pairs_same)
S['class'] = np.select([(S.n_abundant > 0) & (S.n_rare > 0),
                        (S.n_units_detected >= nU / 2) & (S.n_abundant == 0) & (S.n_intermediate == 0),
                        (S.n_abundant > 0) & (S.n_rare == 0) & (S.n_units_detected >= nU / 2)],
                       ['conditionally rare', 'consistently rare', 'consistently abundant/intermediate'], 'other')
S['n_abundant_prev25'] = ((M >= 0.01) & (P >= 0.25)).sum(axis=1)
cr = S['class'] == 'conditionally rare'
S['subtype'] = ''
S.loc[cr & (S.n_rare >= S.n_abundant + S.n_intermediate), 'subtype'] = 'mostly rare, locally abundant'
S.loc[cr & (S.n_rare < S.n_abundant + S.n_intermediate), 'subtype'] = 'mostly abundant, locally rare'
tax = pd.read_csv(tax_f, sep='\t').drop_duplicates('sOTU').set_index('sOTU')
S['top_hit'] = [f"{' '.join(str(tax.loc[s, 'Hit_name']).split()[:2])} ({tax.loc[s, 'Per_ident']:.1f}%)"
                if s in tax.index else 'no BLAST hit' for s in S.index]
yi = pd.read_excel(iso_xlsx, sheet_name='your_isolates')
idm = {'sfa6': 'BF7', 'BF261': 'BF3', 'M1241': 'BF6', 'CSWG105': 'BF5', 'BF170': 'BF2', '78E': 'BF1'}
bf = {idm[r.isolate]: str(r.V4_120bp)[:90] for _, r in yi[yi.V4_copy == 1].iterrows()}
S['BF_isolate'] = S.index.map({v: k for k, v in bf.items()}).fillna('')
print(S['class'].value_counts().to_dict())
print('conditionally rare:', S.loc[cr, 'subtype'].value_counts().to_dict(),
      'prev25-robust', int((cr & (S.n_abundant_prev25 > 0)).sum()),
      'host-contrast', int(S.loc[cr, 'host_contrast_AvsR'].sum()), 'same-host', int(S.loc[cr, 'same_host_contrast_AvsR'].sum()))
for b, s in bf.items():
    print(b, S.loc[s, ['class', 'n_units_detected', 'n_abundant', 'n_intermediate', 'n_rare', 'max_unit_mean_pct']].to_dict() if s in S.index else 'not detected')

# ---- tables
outS = S[['class', 'subtype', 'BF_isolate', 'top_hit', 'n_units_detected', 'n_abundant', 'n_abundant_prev25',
          'n_intermediate', 'n_rare', 'max_unit_mean_pct', 'max_sample_pct', 'median_unit_prev_when_detected',
          'host_contrast_AvsR', 'same_host_contrast_AvsR']].copy()
outS.index.name = 'sOTU_90bp'
outS['abundant_in'] = ['; '.join(C.columns[C.loc[s] == 'A']) for s in outS.index]
outS['rare_in'] = ['; '.join(C.columns[C.loc[s] == 'R']) if outS.loc[s, 'class'] == 'conditionally rare' else '' for s in outS.index]
outS[['class', 'subtype', 'BF_isolate', 'top_hit', 'n_units_detected', 'max_unit_mean_pct']].to_csv(f'{out}/sOTU_classes_all.tsv.gz', sep='\t')
U = info.copy(); U['n_samples'] = u; U.index.name = 'unit'
meth = pd.DataFrame({'note': [
    'Data: amphibian meta-analysis 90-bp V4 deblur sOTUs; samples with >=1,000 reads; skin (external) samples only.',
    'Units: host species x dataset (DOI) x country x life stage with >=15 samples.',
    'Per unit: mean relative abundance (zeros included). Classes: abundant >=1%, intermediate 0.1-1%, rare detected <0.1%, absent = no reads.',
    'Consistently rare: detected in >= half of units and never >=0.1% in any unit.',
    'Conditionally rare: abundant (>=1%) in >=1 unit and rare (<0.1%, detected) in >=1 other unit. Subtype by whether rare units outnumber abundant+intermediate units.',
    'host_contrast_AvsR: abundant vs rare between two host species in the same dataset, country and life stage. same_host_contrast_AvsR: abundant vs rare within one host species across contexts.',
    'n_abundant_prev25: units where abundant AND prevalence >=25% (guards against one-sample blooms).',
    'Pairwise sheet: for each ordered pair of units, % of unit a abundant sOTUs that are abundant/rare/absent in unit b.',
    'top_hit: BLAST top hit (NCBI 16S RefSeq) with % identity over 90 bp; genus labels below ~95% are unreliable.',
    'Caveat: most host species come from one dataset; host and study are separable only within datasets with several co-sampled hosts.']})
cls = outS[outS['class'] != 'other'].sort_values(['class', 'n_abundant'], ascending=[True, False])
with pd.ExcelWriter(f'{out}/conditional_rarity_meta.xlsx') as w:
    meth.to_excel(w, sheet_name='methods', index=False)
    cls.to_excel(w, sheet_name='sOTU_classes')
    outS[outS.BF_isolate != ''].to_excel(w, sheet_name='BF_isolates')
    PW.to_excel(w, sheet_name='pairwise_flips', index=False)
    U.to_excel(w, sheet_name='units')
    M.loc[cls.index].mul(100).round(4).to_excel(w, sheet_name='unit_mean_pct_classified')

# ---- figure
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                     'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 7, 'axes.titlesize': 8})
lab = {'country': 'Same host, country differs', 'stage': 'Same host, life stage differs',
       'country+stage': 'Same host, country + stage differ', 'study+stage': 'Same host, study + stage differ',
       'host only': 'Different host, same dataset/country/stage'}
PW['grp'] = PW.differs.map(lab).fillna('Different host, other factors also differ')
order = list(lab.values()) + ['Different host, other factors also differ']
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.4), gridspec_kw=dict(width_ratios=[1.15, 1], wspace=0.75))
ax = axs[0]; rng = np.random.default_rng(1)
for i, g in enumerate(order):
    v = PW.loc[PW.grp == g, 'rare_or_absent'].values; y = len(order) - 1 - i
    col = '#4477AA' if g.startswith('Same') else '#AA7744'
    ax.scatter(v, y + rng.uniform(-.25, .25, len(v)), s=4 if len(v) > 50 else 14, color=col,
               alpha=.3 if len(v) > 50 else .9, lw=0, rasterized=True)
    ax.plot([np.median(v)] * 2, [y - .33, y + .33], color='k', lw=1.6)
ax.set_yticks(range(len(order)))
ax.set_yticklabels([f'{g}\n(n={int((PW.grp == g).sum())} pairs)' for g in order[::-1]])
ax.set_xlabel('Abundant sOTUs (≥1%) that are rare or\nabsent in the compared group (%)'); ax.set_xlim(-4, 104)
ax.set_title('Abundant taxa become rare more often\nacross sites/stages than across co-sampled hosts', loc='left')
ax.text(-0.02, 1.08, 'a', transform=ax.transAxes, fontsize=10, fontweight='bold', ha='right')
ax = axs[1]
oth = S[S['class'].isin(['other', 'consistently abundant/intermediate'])]
ax.scatter(oth.n_units_detected + rng.uniform(-.3, .3, len(oth)), oth.max_unit_mean_pct, s=2, color='#BBBBBB', lw=0, rasterized=True, label='Other sOTUs')
for c, col in {'consistently rare': '#117777', 'conditionally rare': '#771155'}.items():
    d = S[S['class'] == c]
    ax.scatter(d.n_units_detected + rng.uniform(-.3, .3, len(d)), d.max_unit_mean_pct, s=6, color=col, lw=0, label=f'{c.capitalize()} (n={len(d)})')
ax.set_yscale('log'); ax.axhline(1, color='k', lw=.6, ls='--'); ax.axhline(0.1, color='k', lw=.6, ls=':')
ax.set_ylim(1e-4, 150); ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y, _: f'{y:g}'))
pos = {'BF1': (27, 0.3), 'BF5': (31, 4.5), 'BF3': (38, 0.02), 'BF6': (10, 0.06)}
for b, xy in pos.items():
    s = bf[b]
    if s not in S.index: continue
    x, y = S.loc[s, 'n_units_detected'], S.loc[s, 'max_unit_mean_pct']
    ax.scatter(x, y, marker='*', s=60, color='k', zorder=5)
    ax.annotate(b, (x, y), xytext=xy, fontsize=7, arrowprops=dict(arrowstyle='-', lw=.6), zorder=6, bbox=dict(fc='white', ec='none', pad=.5))
ax.set_xlabel('Host groups where detected (of %d)' % nU); ax.set_ylabel('Highest group mean abundance (%)')
ax.set_title('Most sOTUs stay rare; a few bloom\nin some host groups', loc='left')
ax.legend(loc='upper center', bbox_to_anchor=(0.45, -0.2), ncol=2, frameon=False, markerscale=2, handletextpad=.3, columnspacing=.8)
ax.text(1, 1.2, '1%', fontsize=6); ax.text(1, 0.12, '0.1%', fontsize=6)
ax.text(-0.02, 1.08, 'b', transform=ax.transAxes, fontsize=10, fontweight='bold', ha='right')
fig.savefig(f'{out}/conditional_rarity_meta.png', dpi=300, bbox_inches='tight')
print('done')
