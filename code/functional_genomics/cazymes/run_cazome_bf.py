"""Full dbCAN V13 CAZome + Pfam/NCBIfam panel search on BF1-BF7 Prokka proteomes.
Adapted from the chytrid run_cazome.py / run_hmm.py (same thresholds); CPU capped at 2.
dbCAN thresholds: i-Evalue < 1e-15 if HMM alignment > 80 aa, else < 1e-5; HMM coverage > 0.35
(search run with E=1e-5, domE=1e-5). Pfam models: gathering (GA) cutoffs. NCBIfam/TIGRFAM models: trusted (TC) cutoffs.
Outputs: raw_dbcan_hits.tsv, raw_panel_hits.tsv"""
import pyhmmer, glob, os, time, pandas as pd
from pyhmmer.plan7 import HMMFile
from pyhmmer.easel import SequenceFile, Alphabet
CPU = 2
ISO = [f"BF{i}" for i in range(1, 8)]
abc = Alphabet.amino()
nm = lambda x: x.decode() if isinstance(x, bytes) else x
with HMMFile("db/dbCAN.txt") as f: caz = list(f)
pf, tg = [], []
for p in sorted(glob.glob("db/hmm/*.hmm")):
    with HMMFile(p) as f:
        (pf if os.path.basename(p).startswith("PF") else tg).extend(list(f))
print("HMMs dbCAN", len(caz), "Pfam", len(pf), "NCBIfam/TIGR", len(tg), flush=True)
rc, rp = [], []; t0 = time.time()
for iso in ISO:
    with SequenceFile(f"fg_inputs/{iso}/{iso}.faa", digital=True, alphabet=abc) as sf:
        seqs = sf.read_block()
    plen = {nm(s.name): len(s) for s in seqs}
    for top in pyhmmer.hmmsearch(caz, seqs, cpus=CPU, E=1e-5, domE=1e-5):
        q = top.query; qn = nm(q.name).replace(".hmm", "")
        for hit in top:
            for d in hit.domains:
                a = d.alignment; alen = a.hmm_to - a.hmm_from + 1; cov = alen / q.M
                if (d.i_evalue < 1e-15 if alen > 80 else d.i_evalue < 1e-5) and cov > 0.35:
                    rc.append((iso, nm(hit.name), plen[nm(hit.name)], qn, d.i_evalue, d.score, round(cov, 3), a.hmm_from, a.hmm_to, d.env_from, d.env_to))
    for models, cut, db in ((pf, "gathering", "Pfam"), (tg, "trusted", "NCBIfam")):
        for top in pyhmmer.hmmsearch(models, seqs, cpus=CPU, bit_cutoffs=cut):
            q = top.query
            for hit in top.included:
                for d in hit.domains.included:
                    a = d.alignment
                    rp.append((iso, nm(hit.name), plen[nm(hit.name)], db, nm(q.accession).split(".")[0], nm(q.name), d.i_evalue, d.score, hit.score,
                               round((a.hmm_to - a.hmm_from + 1) / q.M, 3), d.env_from, d.env_to))
    print(iso, len(seqs), round(time.time() - t0), flush=True)
pd.DataFrame(rc, columns=["isolate","protein","prot_len","model","i_evalue","score","hmm_cov","hmm_from","hmm_to","env_from","env_to"]).to_csv("raw_dbcan_hits.tsv", sep="\t", index=False)
pd.DataFrame(rp, columns=["isolate","protein","prot_len","db","accession","model","i_evalue","dom_score","seq_score","hmm_cov","env_from","env_to"]).to_csv("raw_panel_hits.tsv", sep="\t", index=False)
print("done", len(rc), len(rp))
