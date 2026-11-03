"""Single-copy core-genome tree. Usage: 04_core_tree.py <sp> <panaroo_dir> <genomes.tsv> <outdir> [workers]
Per-gene nucleotide MAFFT (--auto, 1 thread each, <=4 concurrent), concatenation (gaps for nothing missing: all genomes present by construction),
FastTree -nt -gtr -gamma."""
import sys, os, subprocess, pandas as pd
from concurrent.futures import ThreadPoolExecutor
sp, pan, gtsv, out = sys.argv[1:5]; W = int(sys.argv[5]) if len(sys.argv) > 5 else 4
genomes = [l.split('\t')[0] for l in open(gtsv) if l.strip()]
core = [l.strip() for l in open(f'{out}/single_copy_core_{sp}.txt') if l.strip()]
pa = pd.read_csv(f'{pan}/gene_presence_absence.csv', low_memory=False).set_index('Gene')
gd = pd.read_csv(f'{pan}/gene_data.csv', usecols=['annotation_id', 'dna_sequence']).set_index('annotation_id').dna_sequence
ad = f'{out}/core_aln_{sp}'; os.makedirs(ad, exist_ok=True)
def run(i_c):
    i, c = i_c
    fa = f'{ad}/g{i:05d}.fa'; al = fa + '.aln'
    if os.path.exists(al) and os.path.getsize(al) > 0: return al
    with open(fa, 'w') as f:
        for g in genomes: f.write(f'>{g}\n{gd[pa.loc[c, g]]}\n')
    with open(al, 'w') as o: subprocess.run(['mafft', '--auto', '--thread', '1', '--quiet', fa], stdout=o, check=True)
    return al
with ThreadPoolExecutor(W) as ex: alns = list(ex.map(run, enumerate(core)))
cat = {g: [] for g in genomes}
for al in alns:
    s = {}; cur = None
    for l in open(al):
        if l.startswith('>'): cur = l[1:].strip(); s[cur] = []
        else: s[cur].append(l.strip().upper())
    for g in genomes: cat[g].append(''.join(s[g]))
L = sum(len(x) for x in cat[genomes[0]])
with open(f'{out}/core_concat_{sp}.fasta', 'w') as f:
    for g in genomes: f.write(f'>{g}\n{"".join(cat[g])}\n')
# variable sites
cols = list(zip(*[''.join(cat[g]) for g in genomes]))
var = sum(1 for c in cols if len(set(x for x in c if x in 'ACGT')) > 1)
with open(f'{out}/core_tree_{sp}.nwk', 'w') as o, open(f'{out}/fasttree_{sp}.log', 'w') as e:
    subprocess.run(['FastTree', '-nt', '-gtr', '-gamma', f'{out}/core_concat_{sp}.fasta'], stdout=o, stderr=e, check=True)
print(f'{sp}: {len(core)} single-copy core genes, alignment {L} bp, {var} variable sites')
