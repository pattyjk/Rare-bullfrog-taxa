"""Post-process a Panaroo pangenome for one species.
Usage: 03_pangenome_analysis.py <species_tag> <isolate_id> <panaroo_dir> <genomes.tsv> <fg_inputs_dir> <outdir>
Outputs: presence/absence matrix (counts), per-isolate specific clusters (<=10% of references),
mapping to Prokka proteins (DIAMOND, >=95% id, >=80% query and subject coverage), islands, core-gene list."""
import sys, os, re, subprocess, pandas as pd, numpy as np
sp, iso, pan, gtsv, fg, out = sys.argv[1:7]
os.makedirs(out, exist_ok=True)
MAXFRAC = 0.10       # isolate-specific: present in <=10% of references
ISL_GAP = 2          # max non-specific genes between specific genes within an island
ISL_MIN = 3          # min specific genes to call a genomic island
genomes = [l.split('\t')[0] for l in open(gtsv) if l.strip()]
refs = [g for g in genomes if g != iso]
pa = pd.read_csv(f'{pan}/gene_presence_absence.csv', low_memory=False)
pa = pa.set_index('Gene')
cnt = pa[genomes].apply(lambda c: c.fillna('').map(lambda x: 0 if x == '' else len(x.split(';'))))
cnt.to_csv(f'{out}/pangenome_presence_absence_{sp}.tsv.gz', sep='\t', compression='gzip')
pres = (cnt > 0).astype(int)
nref = len(refs)
freq = pres[refs].sum(axis=1)
stats = {}
ntot = len(genomes)
k = pres[genomes].sum(axis=1)
stats['n_genomes'] = ntot; stats['n_references'] = nref; stats['total_clusters'] = len(pres)
stats['core_clusters_all_genomes(100%)'] = int((k == ntot).sum())
stats['core_clusters_>=99%'] = int((k >= 0.99 * ntot).sum())
stats['soft_core_95-99%'] = int(((k >= 0.95 * ntot) & (k < 0.99 * ntot)).sum())
stats['shell_15-95%'] = int(((k >= 0.15 * ntot) & (k < 0.95 * ntot)).sum())
stats['cloud_<15%'] = int((k < 0.15 * ntot).sum())
stats['accessory_<99%'] = int((k < 0.99 * ntot).sum())
stats[f'{iso}_clusters'] = int(pres[iso].sum())
stats[f'{iso}_core_>=99%'] = int(((k >= 0.99 * ntot) & (pres[iso] == 1)).sum())
stats[f'{iso}_accessory'] = int(((k < 0.99 * ntot) & (pres[iso] == 1)).sum())
spec = pres.index[(pres[iso] == 1) & (freq <= MAXFRAC * nref)]
stats[f'{iso}_specific_<=10%refs'] = len(spec)
stats[f'{iso}_unique_0refs'] = int(((pres[iso] == 1) & (freq == 0)).sum())
stats['max_refs_allowed_for_specific'] = int(np.floor(MAXFRAC * nref))
# single-copy core
sc = cnt.index[(cnt[genomes] == 1).all(axis=1)]
sc = [c for c in sc if not any(re.search('refound|stop|len|pseudo', str(pa.loc[c, g])) for g in genomes)]
stats['single_copy_core_used_for_tree'] = len(sc)
pd.Series(sc).to_csv(f'{out}/single_copy_core_{sp}.txt', index=False, header=False)
pd.Series(stats).to_csv(f'{out}/pangenome_stats_{sp}.tsv', sep='\t', header=['value'])
# isolate genes: coordinates from the pyrodigal GFF
_r = []
for l in open(f'genes/{iso}.gff'):
    if l.startswith('##FASTA'): break
    if l.startswith('#'): continue
    _r.append(l.rstrip('\n').split('\t'))
gff = pd.DataFrame(_r, columns='ctg src typ start end score strand ph attr'.split())
gff[['start', 'end']] = gff[['start', 'end']].astype(int)
gff = gff[gff.typ == 'CDS'].copy()
gff['lt'] = gff.attr.str.extract(r'ID=([^;]+)')[0]
gff['order'] = range(len(gff))
gff = gff.set_index('lt')
# --- isolate-specific genes: Panaroo (strict) and MMseqs2 (>=50% id, >=80% bidirectional coverage) ---
pan_spec_genes = {}
for c in spec:
    for lt in str(pa.loc[c, iso]).split(';'): pan_spec_genes[lt] = c
mmc = pd.read_csv(f'mm_{sp}_cluster.tsv', sep='\t', names=['rep', 'mem'])
mmc['genome'] = mmc.mem.str.rsplit('_', n=1).str[0]
mP = mmc.groupby(['rep', 'genome']).size().unstack(fill_value=0) > 0
mfreq = mP[refs].sum(axis=1)
mspec = set(mP.index[mP[iso] & (mfreq <= np.floor(MAXFRAC * nref))])
mm_iso = mmc[mmc.genome == iso].set_index('mem')
mm_spec_genes = {m: r for m, r in mm_iso.rep.items() if r in mspec}
stats['mmseqs_total_clusters'] = len(mP); stats['mmseqs_core_100%'] = int(mP[genomes].all(axis=1).sum())
stats[f'mmseqs_{iso}_clusters'] = int(mP[iso].sum()); stats[f'mmseqs_{iso}_specific_<=10%refs'] = len(mspec)
stats[f'mmseqs_{iso}_unique_0refs'] = int((mP[iso] & (mfreq == 0)).sum())
pan_genes = set(';'.join(pa[iso].dropna().astype(str)).split(';'))
allg = sorted(set(pan_spec_genes) | set(mm_spec_genes))
rows = []
for lt in allg:
    pc = pan_spec_genes.get(lt); mc = mm_iso.rep.get(lt)
    rows.append(dict(isolate_locus=lt, cluster=pc if pc else '', panaroo_annotation=pa.loc[pc, 'Annotation'] if pc else '',
        n_refs_with_cluster=int(freq[pc]) if pc else np.nan, frac_refs=round(freq[pc] / nref, 3) if pc else np.nan,
        mmseqs_cluster=mc if mc is not None else '', mmseqs_n_refs=int(mfreq[mc]) if mc is not None else np.nan,
        specific_panaroo=lt in pan_spec_genes, specific_mmseqs=lt in mm_spec_genes, in_panaroo_output=lt in pan_genes))
d = pd.DataFrame(rows)
d['tier'] = np.where(d.specific_panaroo & d.specific_mmseqs, 'both', np.where(d.specific_mmseqs, 'mmseqs_only', 'panaroo_only'))
d['refound'] = d.isolate_locus.str.contains('refound')
d = d.join(gff[['ctg', 'start', 'end', 'strand', 'order']], on='isolate_locus')
d['contig_cov'] = d.ctg.astype(str).str.extract(r'_cov_([0-9.]+)')[0].astype(float)
d['low_coverage_contig_flag'] = d.contig_cov < 5
stats[f'{iso}_specific_union_genes'] = len(d)
stats[f'{iso}_specific_genes_on_lowcov(<5x)_contigs_excluded'] = int(d.low_coverage_contig_flag.sum())
d.loc[d.low_coverage_contig_flag].to_csv(f'{out}/{iso}_specific_lowcov_excluded.tsv', sep='\t', index=False)
d = d[~d.low_coverage_contig_flag].copy()
for t in ['both', 'mmseqs_only', 'panaroo_only']: stats[f'{iso}_specific_genes_tier_{t}'] = int((d.tier == t).sum())
pd.Series(stats).to_csv(f'{out}/pangenome_stats_{sp}.tsv', sep='\t', header=['value'])
# write specific proteins (from pyrodigal faa)
seqs = {}; cur = None
for l in open(f'genes/{iso}.faa'):
    if l.startswith('>'): cur = l[1:].split()[0]; seqs[cur] = ''
    else: seqs[cur] += l.strip()
with open(f'{out}/{iso}_specific.faa', 'w') as f:
    for lt in d.isolate_locus:
        if lt in seqs: f.write(f'>{lt}\n{seqs[lt]}\n')
# map to Prokka proteins
pk = f'{fg}/{iso}/{iso}.faa'
subprocess.run(f'diamond makedb --in {pk} -d {out}/{iso}_prokka --quiet && diamond blastp -q {out}/{iso}_specific.faa -d {out}/{iso}_prokka --threads 4 --quiet --max-target-seqs 1 --id 95 --query-cover 80 --subject-cover 80 --outfmt 6 qseqid sseqid pident qcovhsp scovhsp -o {out}/{iso}_spec_vs_prokka.tsv', shell=True, check=True)
m = pd.read_csv(f'{out}/{iso}_spec_vs_prokka.tsv', sep='\t', names=['isolate_locus', 'prokka_locus_tag', 'prokka_pident', 'qcov', 'scov']).drop_duplicates('isolate_locus')
prod = {}
for l in open(f'{fg}/{iso}/{iso}.gff'):
    if l.startswith('##FASTA'): break
    if l.startswith('#'): continue
    p = l.rstrip('\n').split('\t')
    if len(p) == 9 and p[2] == 'CDS':
        a = dict(x.split('=', 1) for x in p[8].split(';') if '=' in x)
        prod[a.get('locus_tag', a.get('ID'))] = (a.get('product', ''), a.get('gene', ''), a.get('EC_number', ''))
m['prokka_product'] = m.prokka_locus_tag.map(lambda x: prod.get(x, ('', '', ''))[0])
m['prokka_gene'] = m.prokka_locus_tag.map(lambda x: prod.get(x, ('', '', ''))[1])
m['prokka_EC'] = m.prokka_locus_tag.map(lambda x: prod.get(x, ('', '', ''))[2])
d = d.merge(m[['isolate_locus', 'prokka_locus_tag', 'prokka_pident', 'prokka_product', 'prokka_gene', 'prokka_EC']], on='isolate_locus', how='left')
# islands: contiguity on isolate contigs
specset = set(d.isolate_locus)
gl = gff.reset_index().sort_values('order')
isl = {}; iid = 0
for ctg, gg in gl.groupby('ctg', sort=False):
    lts = gg['lt'].tolist(); run = []; last = -99
    for i, lt in enumerate(lts):
        if lt in specset:
            if run and i - last - 1 > ISL_GAP:
                iid += 1; isl.update({x: (iid, len(run)) for x in run}); run = []
            run.append(lt); last = i
    if run: iid += 1; isl.update({x: (iid, len(run)) for x in run})
d['segment'] = d.isolate_locus.map(lambda x: isl.get(x, (None, 0))[0])
d['segment_size'] = d.isolate_locus.map(lambda x: isl.get(x, (None, 0))[1])
seg_ids = {}
for s, n in d.dropna(subset=['segment']).drop_duplicates('segment')[['segment', 'segment_size']].sort_values('segment').values:
    if n >= ISL_MIN: seg_ids[s] = f'{iso}_GI{len(seg_ids)+1:02d}'
d['island_id'] = d.segment.map(lambda s: seg_ids.get(s, ''))
d.loc[d.island_id == '', 'island_id'] = d.loc[d.island_id == '', 'segment'].map(lambda s: f'{iso}_seg{int(s):03d}_small' if pd.notna(s) else 'refound_no_coords')
d = d.drop(columns=['segment'])
d.to_csv(f'{out}/_spec_{iso}.tsv', sep='\t', index=False)
print(pd.Series(stats).to_string())
print('tiers', d.tier.value_counts().to_dict()); print('islands >=3 genes:', len(seg_ids), 'genes in islands:', (d.island_id.str.contains('_GI')).sum(), 'mapped to prokka:', d.prokka_locus_tag.notna().sum(), '/', len(d))
