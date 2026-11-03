"""Island x reference gene-sharing matrix (fraction of island genes whose MMseqs2 50%id/80%cov cluster contains the reference)
and island evidence columns. Usage: 06_island_sharing.py <sp> <iso> <genomes.tsv> <reference_genomes.tsv>"""
import sys, pandas as pd, numpy as np
sp, iso, gtsv, rgt = sys.argv[1:5]
refs = [l.split('\t')[0] for l in open(gtsv) if l.strip() and l.split('\t')[0] != iso]
strain = pd.read_csv(rgt, sep='\t').set_index('accession').strain.fillna('')
acc = lambda x: x.rsplit('_', 1)[0] + '.' + x.rsplit('_', 1)[1]
d = pd.read_csv(f'deliver/isolate_specific_clusters_{iso}.tsv', sep='\t')
c = pd.read_csv(f'mm_{sp}_cluster.tsv', sep='\t', names=['rep', 'mem']); c['genome'] = c.mem.str.rsplit('_', n=1).str[0]
mP = c.groupby(['rep', 'genome']).size().unstack(fill_value=0) > 0
rows = []
for i, sub in d[d.island_id.str.contains('_GI')].groupby('island_id'):
    r = sub.mmseqs_cluster.dropna(); r = r[r.isin(mP.index)]
    rows.append((mP.loc[r, refs].mean(axis=0) if len(r) else pd.Series(0.0, index=refs)).rename(i))
M = pd.DataFrame(rows); M.columns = [acc(x) for x in M.columns]
M.round(3).to_csv(f'deliver/island_gene_sharing_matrix_{iso}.tsv', sep='\t')
isl = pd.read_csv(f'deliver/islands_{iso}.tsv', sep='\t').set_index('island_id')
ev = d.groupby('island_id').agg(n_specific_mmseqs=('specific_mmseqs', 'sum'), n_specific_panaroo=('specific_panaroo', 'sum'))
isl = isl.join(ev)
isl['frac_genes_no_homolog_>=50id_in_>10%refs'] = (isl.n_specific_mmseqs / isl.n_genes).round(2)
isl['most_similar_reference'] = M.idxmax(axis=1).reindex(isl.index)
isl['most_similar_reference_strain'] = isl.most_similar_reference.map(lambda a: strain.get(a, '') if isinstance(a, str) else '')
isl['max_ref_gene_share'] = M.max(axis=1).reindex(isl.index).round(2)
isl['evidence'] = np.where(isl['frac_genes_no_homolog_>=50id_in_>10%refs'] >= 0.5, 'novel_content(>=50% genes lack homologs in >10% refs)', 'mostly_divergent_alleles(homologs >=50% id in refs)')
isl.sort_values('n_genes', ascending=False).to_csv(f'deliver/islands_{iso}.tsv', sep='\t')
print(iso, isl.evidence.value_counts().to_dict())
