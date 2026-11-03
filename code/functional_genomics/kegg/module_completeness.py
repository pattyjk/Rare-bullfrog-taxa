#!/usr/bin/env python3
"""KEGG module completeness from KO sets.
Grammar of KEGG DEFINITION (precedence high->low): '+' / '-' (complex components; '-X' = optional
component), ' ' (AND / sequential steps), ',' (OR / alternatives). Parentheses group.
'--' = step without known KO (excluded from denominator). Embedded module IDs (Mxxxxx) are
evaluated recursively (satisfied if the referenced module is complete).
Completeness = (# top-level steps satisfied) / (# scorable top-level steps). Steps that are entirely
optional (leading '-') are excluded. If the top level itself is an OR of alternative paths,
the best-scoring path is reported. A step is satisfied if any alternative is fully present
(all non-optional complex subunits present).
Multi-line DEFINITIONs are joined with spaces (each line = additional steps).
Usage: module_completeness.py modules_raw_dir ko_presence.tsv(isolate\tKO) out_prefix
"""
import re, sys, os, glob, csv
from collections import defaultdict

def parse_entry(txt):
    d = {'DEFINITION': [], 'CLASS': '', 'NAME': '', 'ENTRY': ''}
    cur = None
    for line in txt.split('\n'):
        if line.startswith('///'): break
        key = line[:12].strip()
        val = line[12:].strip()
        if key: cur = key
        if cur == 'DEFINITION': d['DEFINITION'].append(val)
        elif cur in ('NAME', 'CLASS') and key: d[cur] = val
        elif cur == 'ENTRY' and key: d['ENTRY'] = val.split()[0]
    d['DEFINITION'] = ' '.join(x for x in d['DEFINITION'] if x)
    return d

TOK = re.compile(r'K\d{5}|M\d{5}|--|[(),+\- ]')
def tokenize(s):
    s = re.sub(r'\s+', ' ', s.strip())
    toks = TOK.findall(s)
    if ''.join(toks) != s: raise ValueError('untokenizable: ' + s)
    return toks

class P:
    def __init__(s, toks): s.t = toks; s.i = 0
    def peek(s): return s.t[s.i] if s.i < len(s.t) else None
    def eat(s, x=None):
        tok = s.peek()
        if x is not None and tok != x: raise ValueError(f'expected {x} got {tok} at {s.i}')
        s.i += 1; return tok
    def p_or(s):
        alts = [s.p_and()]
        while s.peek() == ',': s.eat(); alts.append(s.p_and())
        return alts[0] if len(alts) == 1 else ('OR', alts)
    def p_and(s):
        items = [s.p_comp()]
        while s.peek() == ' ':
            s.eat()
            if s.peek() in (None, ')', ','): break
            items.append(s.p_comp())
        return items[0] if len(items) == 1 else ('AND', items)
    def p_comp(s):
        parts = []
        opt = False
        if s.peek() == '-' : s.eat(); opt = True
        parts.append((s.p_fac(), opt))
        while s.peek() in ('+', '-'):
            opt = s.eat() == '-'
            parts.append((s.p_fac(), opt))
        if len(parts) == 1:
            return ('OPT', parts[0][0]) if parts[0][1] else parts[0][0]
        return ('COMP', parts)
    def p_fac(s):
        tok = s.peek()
        if tok == '(':
            s.eat(); n = s.p_or(); s.eat(')'); return n
        if tok == '--': s.eat(); return ('NONE',)
        if tok and (tok[0] in 'KM'): s.eat(); return ('ID', tok)
        raise ValueError(f'bad token {tok}')

def parse(defn):
    p = P(tokenize(defn)); n = p.p_or()
    if p.peek() is not None: raise ValueError('trailing tokens in ' + defn)
    return n

def kos_in(n, acc):
    if n[0] == 'ID': acc.add(n[1])
    elif n[0] in ('OR', 'AND'): [kos_in(c, acc) for c in n[1]]
    elif n[0] == 'COMP': [kos_in(c, acc) for c, _ in n[1]]
    elif n[0] == 'OPT': kos_in(n[1], acc)
    return acc

class Evaluator:
    def __init__(s, trees): s.trees = trees; s.memo = {}
    def ok(s, n, kos):
        t = n[0]
        if t == 'ID':
            x = n[1]
            if x.startswith('K'): return x in kos
            return s.score(x, kos)[2] if x in s.trees else False
        if t == 'OR': return any(s.ok(c, kos) for c in n[1])
        if t == 'AND': return all(s.ok(c, kos) for c in n[1])
        if t == 'COMP': return all(s.ok(c, kos) for c, o in n[1] if not o)
        if t == 'OPT': return True
        if t == 'NONE': return True
    def path_score(s, n, kos):
        steps = n[1] if n[0] == 'AND' else [n]
        pres = tot = 0; missing = []
        for st in steps:
            if st[0] in ('OPT', 'NONE'): continue
            tot += 1
            if s.ok(st, kos): pres += 1
            else: missing.append(fmt(st))
        return pres, tot, missing
    def score(s, mid, kos):
        key = (mid, id(kos))
        if key in s.memo: return s.memo[key]
        s.memo[key] = (0, 1, False, [])  # recursion guard
        tree = s.trees[mid]
        paths = tree[1] if tree[0] == 'OR' else [tree]
        best = max((s.path_score(p, kos) for p in paths), key=lambda r: (r[0] / r[1] if r[1] else 0, r[0]))
        pres, tot, missing = best
        res = (pres, tot, tot > 0 and pres == tot, missing)
        s.memo[key] = res
        return res

def fmt(n):
    t = n[0]
    if t == 'ID': return n[1]
    if t == 'OR': return '(' + ','.join(fmt(c) for c in n[1]) + ')'
    if t == 'AND': return '(' + ' '.join(fmt(c) for c in n[1]) + ')'
    if t == 'COMP': return ''.join(('-' if o else ('+' if i else '')) + fmt(c) for i, (c, o) in enumerate(n[1]))
    if t == 'OPT': return '-' + fmt(n[1])
    return '--'

def main(raw_dir, presence_tsv, out):
    entries = {}
    for f in sorted(glob.glob(os.path.join(raw_dir, 'M*.txt'))):
        e = parse_entry(open(f).read()); entries[e['ENTRY']] = e
    trees, bad = {}, []
    for m, e in entries.items():
        try: trees[m] = parse(e['DEFINITION'])
        except Exception as ex: bad.append((m, str(ex)))
    iso_kos = defaultdict(set)
    with open(presence_tsv) as fh:
        for r in csv.DictReader(fh, delimiter='\t'): iso_kos[r['isolate']].add(r['KO'])
    isolates = sorted(iso_kos)
    ev = Evaluator(trees)
    with open(out + '_long.tsv', 'w') as fo:
        fo.write('isolate\tmodule\tname\tclass\tsteps_present\tsteps_total\tcompleteness\tcomplete\tge075\tKOs_present\tmissing_steps\n')
        wide = {}
        for m in sorted(trees):
            e = entries[m]; allk = kos_in(trees[m], set())
            for iso in isolates:
                pres, tot, comp, miss = ev.score(m, iso_kos[iso])
                frac = pres / tot if tot else 0.0
                wide[(m, iso)] = frac
                kp = sorted(k for k in allk if k in iso_kos[iso] and k.startswith('K'))
                fo.write(f"{iso}\t{m}\t{e['NAME']}\t{e['CLASS']}\t{pres}\t{tot}\t{frac:.3f}\t{int(comp)}\t{int(frac>=0.75)}\t{','.join(kp)}\t{';'.join(miss)}\n")
    with open(out + '_matrix.tsv', 'w') as fo:
        fo.write('module\tname\tclass\t' + '\t'.join(isolates) + '\n')
        for m in sorted(trees):
            e = entries[m]
            fo.write(f"{m}\t{e['NAME']}\t{e['CLASS']}\t" + '\t'.join(f'{wide[(m,i)]:.3f}' for i in isolates) + '\n')
    print(f'modules parsed {len(trees)}; failed {len(bad)}', bad[:5])

if __name__ == '__main__':
    main(*sys.argv[1:4])
