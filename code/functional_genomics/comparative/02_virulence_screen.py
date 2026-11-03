"""fas/att/hyp virulence-locus screen (reference proteins from pFiD188, GenBank JN093097.1, R. fascians D188).
Locus-gene assignment verified by 100% identity of Z29635.1 (fas operon) and AJ311775.1 (att locus) proteins to JN093097 CDS.
Search: DIAMOND blastp (e<=1e-10, all targets) vs pyrodigal proteins of BF7 + all 193 conspecific (ANI>=95%) R. fascians assemblies.
Gene present: hits with >=80% identity whose query coverage summed over ORFs is >=70% (sum handles genes split into >1 ORF, e.g. hyp).
Also: fraction of all 184 pFiD188 CDS present at the same thresholds (pFiD188-like plasmid proxy).
Usage: 02_virulence_screen.py  (after 00 and 01)"""
import io, re, time, subprocess, requests, pandas as pd, numpy as np
from Bio import SeqIO
E_ = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
def fetch(acc):
    t = requests.get(f'{E_}/efetch.fcgi', params={'db': 'nuccore', 'id': acc, 'rettype': 'gbwithparts', 'retmode': 'text'}, timeout=120).text
    open(f'{acc}.gb', 'w').write(t); time.sleep(1); return SeqIO.read(io.StringIO(t), 'genbank')
rec = fetch('JN093097.1')
cds = pd.DataFrame([dict(locus=f.qualifiers.get('locus_tag', [''])[0], prot=f.qualifiers.get('protein_id', [''])[0], aa=f.qualifiers.get('translation', [''])[0]) for f in rec.features if f.type == 'CDS'])
VMAP = {'pFi_052': ('hyp', 'hyp'), 'pFi_064': ('att', 'attR'), 'pFi_065': ('att', 'attX'), 'pFi_066': ('att', 'attA'), 'pFi_067': ('att', 'attB'), 'pFi_068': ('att', 'attC'),
        'pFi_069': ('att', 'attD'), 'pFi_070': ('att', 'attE'), 'pFi_071': ('att', 'attF'), 'pFi_072': ('att', 'attG'), 'pFi_073': ('att', 'attH'), 'pFi_074': ('fas', 'fasR'),
        'pFi_077': ('fas', 'fasA'), 'pFi_078': ('fas', 'fasB'), 'pFi_079': ('fas', 'fasC'), 'pFi_080': ('fas', 'fasD'), 'pFi_81': ('fas', 'fasE'), 'pFi_82': ('fas', 'fasF')}
GL = {'fas': ['fasR', 'fasA', 'fasB', 'fasC', 'fasD', 'fasE', 'fasF'], 'att': ['attR', 'attX', 'attA', 'attB', 'attC', 'attD', 'attE', 'attF', 'attG', 'attH'], 'hyp': ['hyp']}
with open('virulence_ref.faa', 'w') as f:
    for _, r in cds[cds.locus.isin(VMAP)].iterrows(): f.write(f'>{VMAP[r.locus][1]}|{r.locus}|{r.prot}|{VMAP[r.locus][0]}\n{r.aa}\n')
with open('pFiD188_all.faa', 'w') as f:
    for _, r in cds[cds.aa != ''].iterrows(): f.write(f'>{r.locus}\n{r.aa}\n')
with open('orig_loci.faa', 'w') as f:
    for acc in ['Z29635.1', 'AJ311775.1']:
        for ft in fetch(acc).features:
            if ft.type == 'CDS' and 'translation' in ft.qualifiers: f.write(f">{acc}|{ft.qualifiers.get('gene', ['?'])[0]}|{ft.qualifiers['protein_id'][0]}\n{ft.qualifiers['translation'][0]}\n")
rg = pd.read_csv('deliver/reference_genomes.tsv', sep='\t')
rf = rg[(rg.isolate_compared == 'BF7') & rg.status.isin(['used_in_pangenome', 'dereplicated_out(>=99.9%ANI_to_representative)'])]
used = set(rf[rf.status == 'used_in_pangenome'].accession)
ids = ['BF7'] + [a.replace('.', '_') for a in rf.accession]
acc_of = {a.replace('.', '_'): a for a in rf.accession}; acc_of['BF7'] = 'BF7'
sh = lambda c: subprocess.run(c, shell=True, check=True)
sh('cat ' + ' '.join(f'genes/{x}.faa' for x in ids) + ' > rfas_all.faa && diamond makedb --in rfas_all.faa -d rfas_all --threads 4 --quiet')
sh('diamond makedb --in virulence_ref.faa -d vref --quiet && diamond blastp -q orig_loci.faa -d vref --threads 4 --quiet --outfmt 6 qseqid sseqid pident qcovhsp --max-target-seqs 1 -o deliver_check_orig_vs_pFiD188.tsv')
sh('diamond blastp -q virulence_ref.faa -d rfas_all --threads 4 --quiet --max-target-seqs 0 --evalue 1e-10 --outfmt 6 qseqid sseqid pident qcovhsp scovhsp evalue bitscore -o vir_hits.tsv')
sh('diamond blastp -q pFiD188_all.faa -d rfas_all --threads 4 --quiet --max-target-seqs 0 --evalue 1e-10 --outfmt 6 qseqid sseqid pident qcovhsp bitscore -o pfi_hits.tsv')
h = pd.read_csv('vir_hits.tsv', sep='\t', names='q s pid qcov scov ev bits'.split())
h['gene'] = h.q.str.split('|').str[0]; h['genome'] = h.s.str.rsplit('_', n=1).str[0]
agg = h[h.pid >= 80].groupby(['genome', 'gene']).agg(qcov_sum=('qcov', 'sum'), pid_max=('pid', 'max'), hits=('s', lambda x: ';'.join(sorted(x)))).reset_index()
agg['qcov_sum'] = agg.qcov_sum.clip(upper=100)
meta = rg.set_index('accession')[['strain', 'host', 'isolation_source']]
rows = []
for gm in ids:
    sub = agg[agg.genome == gm].set_index('gene')
    for L, gs in GL.items():
        for gn in gs:
            hit = gn in sub.index
            rows.append(dict(genome=acc_of[gm], locus=L, gene=gn, present=int(hit and sub.loc[gn, 'qcov_sum'] >= 70), max_pident=sub.loc[gn, 'pid_max'] if hit else np.nan,
                             summed_qcov=sub.loc[gn, 'qcov_sum'] if hit else np.nan, hit_ORFs=sub.loc[gn, 'hits'] if hit else ''))
vl = pd.DataFrame(rows).merge(meta, left_on='genome', right_index=True, how='left')
vl['in_pangenome_set'] = vl.genome.isin(used) | (vl.genome == 'BF7')
vl['ref_protein'] = vl.gene.map({v[1]: k for k, v in VMAP.items()})
vl = vl[['genome', 'strain', 'host', 'isolation_source', 'in_pangenome_set', 'locus', 'gene', 'ref_protein', 'present', 'max_pident', 'summed_qcov', 'hit_ORFs']]
vl.to_csv('deliver/virulence_loci_presence.tsv', sep='\t', index=False)
p = pd.read_csv('pfi_hits.tsv', sep='\t', names='q s pid qcov bits'.split()); p['genome'] = p.s.str.rsplit('_', n=1).str[0]
p = p[p.pid >= 80].groupby(['genome', 'q']).qcov.sum().reset_index(); p = p[p.qcov >= 70]
pf = p.groupby('genome').q.nunique().reindex(ids, fill_value=0); pf.index = [acc_of[x] for x in pf.index]
W = vl.pivot_table(index='genome', columns='gene', values='present').astype(int)[sum(GL.values(), [])]
for L, gs in GL.items(): W[f'{L}_genes_present(of{len(gs)})'] = W[gs].sum(axis=1)
W['pFiD188_CDS_shared(of184)'] = pf; W['pFiD188_frac'] = (pf / len(cds[cds.aa != ''])).round(3)
W['in_pangenome_set'] = W.index.isin(used) | (W.index == 'BF7')
W['virulence_type'] = np.where(W['fas_genes_present(of7)'] == 7, 'pathogenic-type (complete fas)',
                       np.where((W['att_genes_present(of10)'] > 0) | (W['hyp_genes_present(of1)'] > 0), 'partial (att and/or hyp only, no fas)', 'non-pathogenic-type (no fas/att/hyp)'))
W = W.join(meta); W.to_csv('deliver/virulence_loci_matrix.tsv', sep='\t')
print(W.groupby(['in_pangenome_set', 'virulence_type']).size()); print(W.loc['BF7', ['fas_genes_present(of7)', 'att_genes_present(of10)', 'hyp_genes_present(of1)', 'pFiD188_CDS_shared(of184)']])
