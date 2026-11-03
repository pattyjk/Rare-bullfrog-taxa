"""Fetch NCBI Datasets v2 assembly reports + genomes for R. fascians (taxid 1828; 'Rhodococcus fascians' resolves to the same taxid)
and S. paucimobilis (taxid 13689); QC; skani ANI; greedy 99.9% ANI dereplication.
Usage: 00_fetch_qc_ani_derep.py [--no-download]   (run from the working dir containing fg_inputs/)
Thresholds: exclude NCBI 'contaminated', CheckM completeness <90%, and NCBI 'genome length too small';
species membership ANI >=95% to BF7/BF1 (skani dist); dereplication: greedy centroid, ANI >=99.9% with both align fractions >=50%,
representatives ordered mandatory > assembly level > contig N50; cap 60 per species (not reached)."""
import os, sys, io, json, time, zipfile, shutil, subprocess, requests, pandas as pd, numpy as np
B = 'https://api.ncbi.nlm.nih.gov/datasets/v2'
SP = {'Rfas': ('Rhodococcoides fascians', 'BF7'), 'Spau': ('Sphingomonas paucimobilis', 'BF1')}
MAND = {'Rfas': ['GCF_030813815.1', 'GCF_001894785.1', 'GCF_000761055.1'],   # DSM 43673 (MiGA closest), NBRC 12155 & LMG 3623 (type)
        'Spau': ['GCF_016028435.1', 'GCF_000739895.2', 'GCF_016027095.1']}  # FDAARGOS_881 (MiGA closest), NBRC 13935 & FDAARGOS_908 (type)
LVL = {'Complete Genome': 0, 'Chromosome': 1, 'Scaffold': 2, 'Contig': 3}
def tax(n): return requests.get(f'{B}/taxonomy/taxon/{requests.utils.quote(n)}', timeout=60).json()['taxonomy_nodes'][0]['taxonomy']['tax_id']
def reports(t):
    out, tok = [], None
    while True:
        p = {'page_size': 1000, 'filters.assembly_version': 'current'}
        if tok: p['page_token'] = tok
        j = requests.get(f'{B}/genome/taxon/{t}/dataset_report', params=p, timeout=120).json()
        out += j.get('reports', []); tok = j.get('next_page_token')
        if not tok: return out
assert tax('Rhodococcus fascians') == tax('Rhodococcoides fascians')
rep = {k: reports(tax(v[0])) for k, v in SP.items()}
json.dump(rep, open('ncbi_reports.json', 'w'))
def attr(r, n):
    for a in r.get('assembly_info', {}).get('biosample', {}).get('attributes', []):
        if a['name'] == n: return a['value']
    return ''
rows = []
for sp, v in rep.items():
    for r in v:
        ai = r.get('assembly_info', {}); bs = ai.get('biosample', {}); cm = r.get('checkm_info', {}); st = r.get('assembly_stats', {}); an = r.get('average_nucleotide_identity', {})
        rows.append(dict(species_set=sp, accession=r['accession'], paired=r.get('paired_accession', ''), organism=r['organism']['organism_name'],
            strain=r['organism'].get('infraspecific_names', {}).get('strain', '') or attr(r, 'strain'), culture_collection=attr(r, 'culture_collection'),
            type_material=r.get('type_material', {}).get('type_label', ''), level=ai.get('assembly_level'), refseq_category=ai.get('refseq_category', ''),
            release=ai.get('release_date'), isolation_source=bs.get('isolation_source', ''), host=bs.get('host', ''), geo=bs.get('geo_loc_name', ''),
            atypical=';'.join(ai.get('atypical', {}).get('warnings', [])) if ai.get('atypical') else '',
            checkm_completeness=cm.get('completeness'), checkm_contamination=cm.get('contamination'), ncbi_ani_status=an.get('taxonomy_check_status', ''),
            ncbi_ani_best=an.get('best_ani_match', {}).get('organism_name', ''), size=st.get('total_sequence_length'), contigs=st.get('number_of_contigs'), n50=st.get('contig_n50')))
md = pd.DataFrame(rows)
md['base'] = md.accession.str[4:]; md['pref'] = md.accession.str.startswith('GCF')
md = md.sort_values(['base', 'pref'], ascending=[True, False])
cmx = md.groupby('base').checkm_completeness.max()
u = md.drop_duplicates('base').copy(); u['checkm_completeness'] = u.base.map(cmx)   # one row per GCA/GCF pair, GCF preferred
u['qc_fail'] = ''
u.loc[u.atypical.str.contains('contaminated'), 'qc_fail'] += 'contaminated;'
u.loc[u.checkm_completeness.notna() & (u.checkm_completeness < 90), 'qc_fail'] += 'checkm<90;'
u.loc[u.atypical.str.contains('genome length too small') & (u.qc_fail == ''), 'qc_fail'] += 'genome_length_too_small;'
keep = u[u.qc_fail == '']
os.makedirs('refs', exist_ok=True)
if '--no-download' not in sys.argv:
    todo = [a for a in keep.accession if not os.path.exists(f'refs/{a}.fna')]
    for i in range(0, len(todo), 25):
        r = requests.post(f'{B}/genome/download', json={'accessions': todo[i:i+25], 'include_annotation_type': ['GENOME_FASTA']}, timeout=600); r.raise_for_status()
        z = zipfile.ZipFile(io.BytesIO(r.content))
        for n in z.namelist():
            if n.endswith('.fna'):
                with z.open(n) as f, open(f"refs/{n.split('/')[-2]}.fna", 'wb') as o: shutil.copyfileobj(f, o)
sel, anis, clm = {}, {}, {}
g = keep.set_index('accession')
for sp, (name, iso) in SP.items():
    accs = list(keep[keep.species_set == sp].accession)
    open(f'list_{sp}.txt', 'w').write('\n'.join([f'refs/{a}.fna' for a in accs] + [f'fg_inputs/{iso}/{iso}.fna']) + '\n')
    if not os.path.exists(f'tri_{sp}.tsv'):
        subprocess.run(f'skani triangle -l list_{sp}.txt -t 4 --sparse -o tri_{sp}.tsv', shell=True, check=True)
        subprocess.run(f'skani dist -q fg_inputs/{iso}/{iso}.fna --rl list_{sp}.txt -t 4 -o {iso.lower()}_vs.tsv', shell=True, check=True)
    d = pd.read_csv(f'{iso.lower()}_vs.tsv', sep='\t'); d['acc'] = d.Ref_file.str.extract(r'refs/(.*)\.fna')[0]
    a = d.dropna(subset=['acc']).set_index('acc').reindex(accs); anis[sp] = a
    members = list(a.index[a.ANI >= 95])
    t = pd.read_csv(f'tri_{sp}.tsv', sep='\t')
    t['r'] = t.Ref_file.str.extract(r'refs/(.*)\.fna')[0]; t['q'] = t.Query_file.str.extract(r'refs/(.*)\.fna')[0]
    t = t.dropna(subset=['r', 'q']); t = t[t.r.isin(members) & t.q.isin(members)]
    S = {}
    for r_, q_, v, f1, f2 in zip(t.r, t.q, t.ANI, t.Align_fraction_ref, t.Align_fraction_query):
        if f1 >= 50 and f2 >= 50: S[(r_, q_)] = S[(q_, r_)] = v
    order = sorted(members, key=lambda x: (x not in MAND[sp], LVL.get(g.loc[x, 'level'], 4), -g.loc[x, 'n50']))
    reps, cl = [], {}
    for x in order:
        hit = [r for r in reps if S.get((x, r), 0) >= 99.9]
        if hit: cl[x] = hit[0]
        else: reps.append(x); cl[x] = x
    assert len(reps) <= 60
    for m in MAND[sp]:
        if m in members and m not in reps: reps.append(m)
    sel[sp], clm[sp] = reps, cl
    with open(f'genomes_{sp}.tsv', 'w') as f:
        f.write(f'{iso}\tfg_inputs/{iso}/{iso}.fna\n')
        for a_ in reps: f.write(f"{a_.replace('.', '_')}\trefs/{a_}.fna\n")
    print(sp, 'assemblies', (u.species_set == sp).sum(), 'pass QC', len(accs), 'ANI>=95', len(members), 'used', len(reps))
def status(r):
    if r.qc_fail: return 'excluded_QC:' + r.qc_fail.rstrip(';')
    a = anis[r.species_set].ANI.get(r.accession, np.nan)
    if pd.isna(a): return 'excluded_ANI_not_reported(<~80%)'
    if a < 95: return 'excluded_ANI<95'
    return 'used_in_pangenome' if r.accession in sel[r.species_set] else 'dereplicated_out(>=99.9%ANI_to_representative)'
u['status'] = u.apply(status, axis=1)
for col, src in [('ANI_to_isolate', 'ANI'), ('AF_ref', 'Align_fraction_ref'), ('AF_isolate', 'Align_fraction_query')]:
    u[col] = u.apply(lambda r: anis[r.species_set][src].get(r.accession, np.nan) if r.qc_fail == '' else np.nan, axis=1)
u['derep_representative'] = u.apply(lambda r: clm[r.species_set].get(r.accession, ''), axis=1)
u['mandatory_keep'] = u.apply(lambda r: r.accession in MAND[r.species_set], axis=1)
u['isolate_compared'] = u.species_set.map({k: v[1] for k, v in SP.items()})
u['species_target'] = u.species_set.map({'Rfas': 'Rhodococcoides fascians (BF7 comparison)', 'Spau': 'Sphingomonas paucimobilis (BF1 comparison)'})
cols = ['isolate_compared', 'species_target', 'accession', 'paired', 'organism', 'strain', 'culture_collection', 'type_material', 'mandatory_keep', 'level', 'refseq_category', 'release', 'size', 'contigs', 'n50',
        'checkm_completeness', 'checkm_contamination', 'atypical', 'ncbi_ani_status', 'ncbi_ani_best', 'isolation_source', 'host', 'geo', 'ANI_to_isolate', 'AF_ref', 'AF_isolate', 'derep_representative', 'status']
os.makedirs('deliver', exist_ok=True)
u[cols].sort_values(['isolate_compared', 'status', 'ANI_to_isolate'], ascending=[True, False, False]).to_csv('deliver/reference_genomes.tsv', sep='\t', index=False)
