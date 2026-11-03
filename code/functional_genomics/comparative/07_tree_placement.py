"""Isolate placement: core-SNP distances (pairwise, ACGT-only sites) from the concatenated single-copy core alignment,
sister clade on the midpoint-rooted FastTree tree, labelled Newick (accession|strain|fas status for R. fascians).
Usage: 07_tree_placement.py <sp> <iso>"""
import sys, numpy as np, pandas as pd
from Bio import SeqIO, Phylo
sp, iso = sys.argv[1:3]
rg = pd.read_csv('deliver/reference_genomes.tsv', sep='\t').set_index('accession')
lab = lambda x: x if x == iso else f"{x.rsplit('_', 1)[0]}.{x.rsplit('_', 1)[1]}"
aln = {r.id: np.frombuffer(str(r.seq).encode(), dtype='S1') for r in SeqIO.parse(f'out_{sp}/core_concat_{sp}.fasta', 'fasta')}
ids = list(aln); A = np.vstack([aln[i] for i in ids]); ok = np.isin(A, [b'A', b'C', b'G', b'T'])
D = pd.DataFrame(0, index=[lab(i) for i in ids], columns=[lab(i) for i in ids])
for i in range(len(ids)):
    for j in range(i + 1, len(ids)):
        m = ok[i] & ok[j]; D.iloc[i, j] = D.iloc[j, i] = int((A[i][m] != A[j][m]).sum())
D.to_csv(f'deliver/core_snp_matrix_{iso}_set.tsv', sep='\t')
D[iso].drop(iso).sort_values().rename(f'core_SNPs_to_{iso}').to_csv(f'deliver/core_snp_distance_to_{iso}.tsv', sep='\t', header=True)
vt = {}
if sp == 'Rfas':
    W = pd.read_csv('deliver/virulence_loci_matrix.tsv', sep='\t', index_col=0)
    vt = W.virulence_type.map({'pathogenic-type (complete fas)': 'fas+', 'partial (att and/or hyp only, no fas)': 'att/hyp-only', 'non-pathogenic-type (no fas/att/hyp)': 'fas-'}).to_dict()
t = Phylo.read(f'out_{sp}/core_tree_{sp}.nwk', 'newick'); t.root_at_midpoint()
leaf = [c for c in t.get_terminals() if c.name == iso][0]; par = t.get_path(leaf)[-2]
print(iso, 'sister:', [lab(c.name) for c in par.get_terminals() if c.name != iso], 'nearest:', D[iso].drop(iso).sort_values().head(3).to_dict(), 'aln bp', A.shape[1])
t0 = Phylo.read(f'out_{sp}/core_tree_{sp}.nwk', 'newick')
for c in t0.get_terminals():
    a = lab(c.name); st = rg.strain.get(a, '') if a in rg.index else ''
    st = '' if pd.isna(st) else str(st)
    nm = a + ('|' + st.replace(' ', '_') if st else '') + (('|' + vt.get(a, '')) if vt else '')
    c.name = nm.replace(':', '_').replace(',', '_').replace('(', '').replace(')', '')
Phylo.write(t0, f'deliver/core_tree_{iso}_{sp}_labeled.nwk', 'newick')
Phylo.write(Phylo.read(f'out_{sp}/core_tree_{sp}.nwk', 'newick'), f'deliver/core_tree_{iso}_{sp}.nwk', 'newick')
