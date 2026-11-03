#!/usr/bin/env python3
"""Top families/genera of consistently and conditionally rare sOTUs vs all other sOTUs.

Inputs : sOTU_classes_all.tsv.gz (from conditional_rarity.py), meta_analysis_taxonomy.txt
Outputs: rare_classes_taxonomy.png, rare_classes_taxonomy.xlsx
Rank rule: BLAST top-hit lineage; family if identity >=86.5%, genus if >=94.5% (Yarza et al. 2014), else Unclassified.
Usage  : python rare_classes_taxonomy.py sOTU_classes_all.tsv.gz TAX.txt OUTDIR
"""
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
cls_f, tax_f, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
S = pd.read_csv(cls_f, sep='\t', index_col=0).fillna({'subtype': '', 'BF_isolate': ''})
tax = pd.read_csv(tax_f, sep='\t').drop_duplicates('sOTU').set_index('sOTU')
def ranks(s):
    if s not in tax.index: return ('Unclassified', 'Unclassified')
    r = tax.loc[s]; toks = [x.strip() for x in str(r.Taxonomy if pd.notna(r.Taxonomy) else '').split(';') if x.strip()]
    fam = next((x for x in toks if x.endswith('aceae')), None)
    g = toks[-1] if toks and toks[-1][0].isupper() and ' ' not in toks[-1] and not toks[-1].endswith(('aceae', 'ales', 'incertae sedis')) else None
    return (fam if fam and r.Per_ident >= 86.5 else 'Unclassified', g if g and r.Per_ident >= 94.5 else 'Unclassified')
S[['family', 'genus']] = pd.DataFrame([ranks(s) for s in S.index], index=S.index)
grp = {'Consistently rare': S['class'] == 'consistently rare', 'Conditionally rare': S['class'] == 'conditionally rare',
       'All other sOTUs': S['class'].isin(['other', 'consistently abundant/intermediate'])}
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'axes.titlesize': 8})
colc = {'Consistently rare': '#117777', 'Conditionally rare': '#771155'}
bg = S[grp['All other sOTUs']]
fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.4), gridspec_kw=dict(hspace=0.42, wspace=0.85))
res = {}
for i, lvl in enumerate(['family', 'genus']):
    bgp = bg[lvl].value_counts(normalize=True) * 100
    plural = {'family': 'families', 'genus': 'genera'}[lvl]
    for j, (g, col) in enumerate(colc.items()):
        d = S[grp[g]]; n = len(d); top = d.loc[d[lvl] != 'Unclassified', lvl].value_counts().head(15)
        pct = top / n * 100; ax = axs[i, j]; y = np.arange(len(top))[::-1]
        ax.barh(y, pct.values, color=col, height=.7, label=f'{g} sOTUs')
        ax.scatter(bgp.reindex(top.index).fillna(0).values, y, marker='|', s=60, color='k', lw=1.2, zorder=3, label='All other sOTUs')
        ax.set_yticks(y); ax.set_yticklabels(top.index, style='italic' if lvl == 'genus' else 'normal')
        for yy, (nm, c) in zip(y, top.items()): ax.text(pct[nm] + 0.15, yy, str(c), va='center', fontsize=6)
        ax.set_xlim(0, max(pct.max(), bgp.reindex(top.index).max()) * 1.2)
        ax.set_title(f'{g}: top {plural}\n(n={n}; {(d[lvl] == "Unclassified").mean() * 100:.0f}% unclassified)', loc='left')
        if i == 1: ax.set_xlabel('% of sOTUs in class')
        ax.text(-0.02, 1.12, 'abcd'[i * 2 + j], transform=ax.transAxes, fontsize=10, fontweight='bold', ha='right')
        res[(lvl, g)] = pd.DataFrame({'n_sOTUs': top, 'pct_of_class': pct.round(2), 'pct_all_other_sOTUs': bgp.reindex(top.index).fillna(0).round(2)})
        print(lvl, g, n, top.head(5).to_dict())
axs[0, 1].legend(loc='lower right', frameon=False, fontsize=6)
fig.savefig(f'{out}/rare_classes_taxonomy.png', dpi=300, bbox_inches='tight')
with pd.ExcelWriter(f'{out}/rare_classes_taxonomy.xlsx') as w:
    pd.DataFrame({'note': ['Taxonomy = BLAST top hit (NCBI 16S RefSeq) of each 90-bp sOTU.',
                           'Family if identity >=86.5%, genus if >=94.5% (Yarza et al. 2014), else Unclassified.',
                           'Classes from conditional_rarity.py.',
                           'pct_all_other_sOTUs = share among sOTUs in neither class (background).',
                           'Many genera tie at 2 sOTUs among conditionally rare; the top-15 cut among ties is arbitrary.']}).to_excel(w, sheet_name='methods', index=False)
    for (lvl, g), d in res.items(): d.to_excel(w, sheet_name=f'{g.split()[0]}_{lvl}')
    for g in ['Consistently rare', 'Conditionally rare']:
        S.loc[grp[g], ['family', 'genus', 'top_hit', 'subtype', 'BF_isolate']].to_excel(w, sheet_name=f'{g.split()[0]}_all_sOTUs')
print('done')
