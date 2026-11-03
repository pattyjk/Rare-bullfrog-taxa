"""Consistent gene calling (pyrodigal, single-genome training, table 11) on all genomes.
Writes Prokka-style GFF3 (+##FASTA) for Panaroo, plus .faa/.ffn. Usage: 01_call_genes.py genomes.tsv outdir  (tsv: genome_id<TAB>fasta)"""
import sys, os, pyrodigal
from multiprocessing import Pool
from Bio import SeqIO
def run(args):
    gid, fa, out = args
    if os.path.exists(f'{out}/{gid}.gff'): return gid
    recs=[(r.id,str(r.seq).upper()) for r in SeqIO.parse(fa,'fasta')]
    gf=pyrodigal.GeneFinder(meta=False)
    gf.train(*[s.encode() for _,s in recs], translation_table=11)
    gff=open(f'{out}/{gid}.gff','w'); faa=open(f'{out}/{gid}.faa','w'); ffn=open(f'{out}/{gid}.ffn','w')
    gff.write('##gff-version 3\n')
    for cid,s in recs: gff.write(f'##sequence-region {cid} 1 {len(s)}\n')
    n=0
    for cid,s in recs:
        for gn in gf.find_genes(s.encode()):
            n+=1; lt=f'{gid}_{n:05d}'
            strand='+' if gn.strand==1 else '-'
            gff.write(f'{cid}\tpyrodigal\tCDS\t{gn.begin}\t{gn.end}\t.\t{strand}\t0\tID={lt};locus_tag={lt};product=hypothetical protein\n')
            prot=gn.translate(translation_table=11).rstrip('*')
            faa.write(f'>{lt} {cid}\n{prot}\n'); ffn.write(f'>{lt} {cid}\n{gn.sequence()}\n')
    gff.write('##FASTA\n')
    for cid,s in recs:
        gff.write(f'>{cid}\n'+'\n'.join(s[i:i+60] for i in range(0,len(s),60))+'\n')
    for f in (gff,faa,ffn): f.close()
    return gid
if __name__=='__main__':
    tsv,out=sys.argv[1],sys.argv[2]; os.makedirs(out,exist_ok=True)
    jobs=[(l.split('\t')[0],l.split('\t')[1].strip(),out) for l in open(tsv) if l.strip()]
    with Pool(int(sys.argv[3]) if len(sys.argv)>3 else 4) as p:
        for g in p.imap_unordered(run,jobs): print(g,flush=True)
