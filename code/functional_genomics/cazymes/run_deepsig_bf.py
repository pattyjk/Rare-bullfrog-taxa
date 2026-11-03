"""Run DeepSig (Savojardo et al. 2018; GitHub BolognaBiocomp/deepsig commit 727941b) signal-peptide
detection on macOS. The bundled CRF cleavage binary (tools/biocrf-static) is Linux-only, so the DeepSig
DNN ensemble (20 models, unchanged thresholds: gramn 0.7575, gramp 0.72, euk 0.725) is called directly and
cleavage sites are decoded with DeepSig's own pure-Python CRF (crfdecoding.py, same models/decoding).
Usage: python run_deepsig_bf.py <fasta> <euk|gramp|gramn> <out.tsv>"""
# NOTE (2026-11-02): in this macOS environment the pure-Python CRF fallback decoded 0 of 20 models, so NO cleavage
# sites were produced (cleavage_pos = '-' for all 31,109 proteins). Only the DNN signal-peptide/transmembrane/other
# class is reported. See outputs/08_functional_genomics/methods/methods_cazymes.md.

import os, sys
os.environ["DEEPSIG_ROOT"] = os.path.abspath("deepsig"); sys.path.insert(0, os.environ["DEEPSIG_ROOT"])
os.environ["OMP_NUM_THREADS"] = "2"; os.environ["TF_NUM_INTEROP_THREADS"] = "1"; os.environ["TF_NUM_INTRAOP_THREADS"] = "2"
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
import numpy, tensorflow as tf
tf.config.threading.set_intra_op_parallelism_threads(2); tf.config.threading.set_inter_op_parallelism_threads(1)
import deepsig.deepsigconfig as cfg
from deepsig.helpers import readdata, detectsp, runCRF2
fa, org, out = sys.argv[1:4]
X, accs, seqs = readdata(fa, cfg.NTERM)
Y, Ytm, Ync, cls, Ytmn, Yncn = detectsp(X, org)
idx = [i for i in range(len(cls)) if cls[i] == 2]
cleav = {i: "-" for i in range(len(cls))}
if idx:
    P = numpy.array([numpy.concatenate((X[i], numpy.zeros((cfg.NTERM, 1))), axis=1) for i in idx])
    C = []
    for i in range(5):
        for j in range(4):
            m = os.path.join(cfg.DEEPSIG_ROOT, cfg.CRF_MODEL_DIR, org, "model.seq.w%d.s%s.l%d.%d.%d" % (cfg.CRF_WINDOWS[org], cfg.CRF_PARAMS[org][i]['sigma'], cfg.NTERM, i, j))
            try:
                pred, _ = runCRF2(P, m, cfg.CRF_WINDOWS[org], cfg.CRF_PARAMS[org][i]['decoding'])
                C.append(pred)
            except Exception:   # pure-Python posterior-Viterbi decoding fails numerically; skip that model
                pass
    print("CRF models decoded:", len(C), "of 20")
    if C:
        C = numpy.mean(numpy.array(C), axis=0)
        for k, i in enumerate(idx): cleav[i] = int((C[k] >= 0.5).sum())
names = {2: "SignalPeptide", 1: "Transmembrane", 0: "Other"}
with open(out, "w") as f:
    f.write("protein\tdeepsig_class\tsp_prob\treliability\tcleavage_pos\n")
    for i, a in enumerate(accs):
        rel = Y[i] if cls[i] == 2 else (Ytmn[i] if cls[i] == 1 else Yncn[i])
        f.write(f"{a}\t{names[cls[i]]}\t{Y[i]:.3f}\t{rel:.3f}\t{cleav[i]}\n")
print("done", len(accs), "SP", len(idx))
