"""Pfam annotation (hmmsearch --cut_ga, Pfam-A 38.2) of isolate-specific proteins; contig coverage flags; BGC-domain flags per island.
Usage: 05_pfam_islands.py <iso> <outdir> <fg_inputs_dir>"""
import sys, subprocess, pandas as pd, re
iso, out, fg = sys.argv[1:4]
subprocess.run(f'hmmsearch --cut_ga --cpu 3 --noali --domtblout {out}/{iso}_specific_pfam.domtbl pfam/Pfam-A.hmm {out}/{iso}_specific.faa > /dev/null', shell=True, check=True)
rows = []
for l in open(f'{out}/{iso}_specific_pfam.domtbl'):
    if l.startswith('#'): continue
    p = l.split()
    rows.append(dict(isolate_locus=p[0], pfam=p[3], pfam_acc=p[4].split('.')[0], ievalue=float(p[12]), env_from=int(p[19])))
pf = pd.DataFrame(rows, columns=['isolate_locus','pfam','pfam_acc','ievalue','env_from'])
pf = pf.sort_values(['isolate_locus', 'env_from'])
pfs = pf.groupby('isolate_locus').agg(pfam_domains=('pfam', lambda x: ';'.join(x)), pfam_accessions=('pfam_acc', lambda x: ';'.join(dict.fromkeys(x))))
d = pd.read_csv(f'{out}/_spec_{iso}.tsv', sep='\t').merge(pfs, left_on='isolate_locus', right_index=True, how='left')
# contig coverage from SPAdes header
d['contig_len'] = d.ctg.astype(str).str.extract(r'_length_([0-9]+)_')[0].astype(float)
# BGC-indicative domains (antiSMASH core-domain proxies)
BGC = {'NRPS': ['Condensation', 'AMP-binding', 'PP-binding', 'NRPS'], 'PKS_typeI/II': ['ketoacyl-synt', 'Ketoacyl-synt_C', 'KAsynt_C_assoc', 'Chal_sti_synt_N', 'Chal_sti_synt_C'],  # KR/DH/AT alone are not used (also hit non-PKS enzymes)
       'terpene': ['Terpene_synth', 'Terpene_synth_C', 'SQHop_cyclase_N', 'SQHop_cyclase_C', 'Lycopene_cycl', 'TPS'], 'lanthipeptide/RiPP': ['LANC_like', 'Lant_dehydr_N', 'Lant_dehydr_C', 'DUF4135', 'Pkinase_fungal_LanL'],
       'NRPS-independent_siderophore': ['IucA_IucC', 'FhuF'], 'betalactone': ['DUF1999'], 'ectoine': ['Ectoine_synth'], 'butyrolactone': ['AfsA'], 'arylpolyene/other': ['Thiopep_pre', 'YcaO']}
def bgc(doms):
    if pd.isna(doms): return ''
    ds = set(doms.split(';')); hits = [k for k, v in BGC.items() if ds & set(v)]
    return ';'.join(hits)
d['bgc_domain_class'] = d.pfam_domains.map(bgc)
isl = d.groupby('island_id').agg(n=('isolate_locus', 'size'), bgc=('bgc_domain_class', lambda x: ';'.join(sorted(set(';'.join(x).split(';')) - {''}))),
                                 ampb=('pfam_domains', lambda x: sum('AMP-binding' in str(v) for v in x)), cond=('pfam_domains', lambda x: sum('Condensation' in str(v) for v in x)),
                                 ks=('pfam_domains', lambda x: sum('ketoacyl-synt' in str(v) for v in x)))
isl['pcp'] = d.groupby('island_id').pfam_domains.agg(lambda x: sum('PP-binding' in str(v) for v in x))
def flag(r):
    if r.cond >= 1 and r.ampb >= 1: return 'likely_NRPS_BGC(C+A domains)'
    if r.ks >= 1 and r.ampb >= 1 and r.pcp >= 1: return 'likely_PKS/NRPS_hybrid_BGC(KS+A+PCP)'
    if r.ampb >= 1 and r.pcp >= 1: return 'possible_NRPS-like(A+PCP, no C domain)'
    if r.ks >= 1 and 'PKS' in r.bgc: return 'likely_PKS_BGC(KS domain)'
    if r.bgc: return 'BGC-associated_domains:' + r.bgc
    return ''
isl['island_bgc_flag'] = isl.apply(flag, axis=1)
import numpy as np
_c = []
for l in open(f'{fg}/{iso}/{iso}.fna'):
    if l.startswith('>'):
        m = re.search(r'_length_([0-9]+)_cov_([0-9.]+)', l)
        if m: _c.append((int(m.group(1)), float(m.group(2))))
_c = pd.DataFrame(_c, columns=['len', 'cov']); _c = _c[_c.len >= 10000].sort_values('cov')
chrom_cov = float(_c['cov'].values[np.searchsorted(_c.len.cumsum().values, _c.len.sum() / 2)])
d['cov_rel_to_median'] = (d.contig_cov / chrom_cov).round(2)
d['putative_plasmid_contig(cov>=1.5x_median)'] = d.cov_rel_to_median >= 1.5
print('length-weighted median coverage of contigs >=10kb:', round(chrom_cov, 2))
d = d.merge(isl[['island_bgc_flag']], left_on='island_id', right_index=True, how='left')
d.loc[~d.island_id.str.contains('_GI'), 'island_bgc_flag'] = d.loc[~d.island_id.str.contains('_GI'), 'bgc_domain_class'].map(lambda x: ('gene-level:' + x) if x else '')
cols = ['island_id', 'segment_size', 'cluster', 'isolate_locus', 'ctg', 'start', 'end', 'strand', 'contig_len', 'contig_cov', 'cov_rel_to_median', 'putative_plasmid_contig(cov>=1.5x_median)', 'low_coverage_contig_flag', 'tier', 'specific_panaroo', 'specific_mmseqs', 'in_panaroo_output', 'n_refs_with_cluster', 'frac_refs', 'mmseqs_cluster', 'mmseqs_n_refs',
        'prokka_locus_tag', 'prokka_pident', 'prokka_gene', 'prokka_product', 'prokka_EC', 'panaroo_annotation', 'pfam_domains', 'pfam_accessions', 'bgc_domain_class', 'island_bgc_flag', 'refound']
d = d.sort_values(['ctg', 'start'])[cols].rename(columns={'ctg': 'contig', 'segment_size': 'n_specific_genes_in_segment', 'isolate_locus': 'pyrodigal_locus'})
d.to_csv(f'deliver/isolate_specific_clusters_{iso}.tsv', sep='\t', index=False)
# island summary
s = d[d.island_id.str.contains('_GI')].groupby('island_id').agg(contig=('contig', 'first'), start=('start', 'min'), end=('end', 'max'), n_genes=('pyrodigal_locus', 'size'),
     n_tier_both=('tier', lambda x: int((x == 'both').sum())), n_mmseqs_0refs=('mmseqs_n_refs', lambda x: int((x == 0).sum())), contig_cov=('contig_cov', 'first'), cov_rel_to_median=('cov_rel_to_median', 'first'), putative_plasmid=('putative_plasmid_contig(cov>=1.5x_median)', 'first'),
     bgc_flag=('island_bgc_flag', 'first'), products=('prokka_product', lambda x: '; '.join(pd.Series([v for v in x if isinstance(v, str) and v != 'hypothetical protein']).value_counts().index[:8])),
     n_hypothetical=('prokka_product', lambda x: int(sum((not isinstance(v, str)) or v == 'hypothetical protein' for v in x))))
s['length_bp'] = s.end - s.start + 1
s.to_csv(f'deliver/islands_{iso}.tsv', sep='\t')
print(iso, 'specific genes', len(d), 'pfam-annotated', d.pfam_domains.notna().sum(), 'islands', len(s), 'low-cov genes', int(d.low_coverage_contig_flag.sum()))
print(s.sort_values('n_genes', ascending=False).head(12)[['contig', 'n_genes', 'length_bp', 'contig_cov', 'bgc_flag', 'products']].to_string())
