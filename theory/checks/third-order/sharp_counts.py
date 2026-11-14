# Part IV, Section 16 (Lemma mg:lem:smallcounts, Corollary mg:cor:sharp): the exact counts that sharpen the constants, each made by two programs.
# Usage: python3 sharp_counts.py      (about 10 minutes; writes orders/sharp_counts.json)
"""(i) Rooted sets: S_n(root) is the number of sets P of n plaquettes of Z^D (D = d-1) with every plaquette joined to the vertices of the root front by a
chain of plaquettes of P sharing vertices (the sets of the rooted counting in the proofs of Lemmas mg:lem:twofrontcount and mg:lem:jsums6). The root is a
plaquette boundary (four edges) or a six-edge front (Z^2: the flat domino; Z^3: the flat and the bent domino and the skew hexagon, the largest count used).
Program 1: rooted_sets.c (Redelmeier's method); program 2: rooted_rs.c (reverse search: the parent of P is P minus its largest removable plaquette).
Sizes: n <= 6 (four-edge root) and n <= 5 (six-edge) in Z^2, n <= 5 and n <= 4 in Z^3. Such sets span at most 7 lattice units in every direction, so
the counts hold on the tori with sides >= 8 as in infinite volume.
(ii) One-plaquette junctions with three or more front-pairs (Lemma mg:lem:jsums3): the candidate junctions with one plaquette containing p = (bd q, s-1)
and further fronts of the sizes K. The classes are those of the tree bound with a-exponent E <= 5 for p' of four edges, (4; 8), (4; 9), (4; 10),
(4; 4,4), (4; 4,6), and, as checks of Lemma mg:lem:betafour(iii) and its extension, (4,4), (4,6), (4,7), (6,6), (4,4,7), (4,4,8), (4,6,6).
Program 1: order_counts.count_junctions (Python); program 2: junctions.c (C).
(iii) Closed polymers through a plaquette of Z^3 (d=3) with 10 <= n <= 16 plaquettes (Lemma mg:lem:theta7 counts them to 14):
closed_sets.c (C) and closed_sets_ext.py (Python)."""
import os, sys, json, subprocess, tempfile, re
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from order_counts import count_junctions

tmp = tempfile.mkdtemp()
def build(src):
    out = os.path.join(tmp, os.path.splitext(src)[0])
    subprocess.run(['cc', '-O2', '-o', out, os.path.join(here, src), '-lm'], check=True)
    return out
RS, RR, JC = build('rooted_sets.c'), build('rooted_rs.c'), build('junctions.c')
def run(cmd):
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
def counts_from(line):
    return [int(x) for x in line.strip().split(':')[-1].split()]

ok = True
rooted = {}
SIZES = {2: {'q': 6, 'flat': 5}, 3: {'q': 5, 'flat': 4, 'bent': 4, 'skew': 4}}
for D, shapes in SIZES.items():
    d = D + 1; six = {}
    for shape, nmax in shapes.items():
        a = counts_from(run([RS, str(D), str(nmax), shape]).splitlines()[-1])
        b = counts_from(run([RR, str(D), str(nmax), shape]).splitlines()[-1])
        same = a == b; ok = ok and same
        print(f"Z^{D}, root {shape}: S_1..S_{nmax} = {a} (Redelmeier) == {b} (reverse search): {same}", flush=True)
        if shape == 'q': rooted.setdefault(str(d), {})['4'] = {str(n + 1): v for n, v in enumerate(a)}
        else: six[shape] = a
    L = min(len(v) for v in six.values())
    rooted[str(d)]['6'] = {str(n + 1): max(v[n] for v in six.values()) for n in range(L)}
    print(f"Z^{D}: six-edge root, the largest over the shapes: {rooted[str(d)]['6']}", flush=True)

CLASSES = [(4, 4), (4, 6), (4, 7), (6, 6), (4, 8), (4, 9), (4, 10), (4, 4, 4), (4, 4, 6), (4, 4, 7), (4, 4, 8), (4, 6, 6)]
oneplaq = {}
for D in (2, 3):
    d = D + 1; terms = []; empty = []
    for K in CLASSES:
        a = dict(count_junctions(D, 1, K)).get(1, 0)
        out = run([JC, str(D), '1', ','.join(map(str, K))])
        m = re.search(r'junctions by n:(.*)', out)
        got = dict(tuple(int(y) for y in x.split(':')) for x in m.group(1).split()) if (m and m.group(1).strip()) else {}
        b = got.get(1, 0)
        same = a == b; ok = ok and same
        print(f"d={d}, one plaquette, further fronts {K}: {a} (order_counts) == {b} (junctions.c): {same}", flush=True)
        if a == 0: empty.append(list(K))
        elif K[0] == 4 and len(K) == 2:              # the pairs with p' the four-edge front: E = K[1]/2 (2E = K[1]) and one factor e^{c_J}
            terms.append([a, K[1], 4, 1])
        else: raise SystemExit(f'unexpected non-empty class {K}')
    oneplaq[str(d)] = {'E0': 5, 'terms': terms, 'empty': empty}
# (iii) closed polymers through a plaquette of Z^3 (d=3), n <= CMAX: closed_sets.c (C) and closed_sets_ext.py (Python), the counts of Lemma mg:lem:theta7
CMAX = 16
CS = build('closed_sets.c')
a = run([CS, '3', str(CMAX)]).strip().splitlines()[-1]
ca = {int(x.split(':')[0]): int(x.split(':')[1]) for x in a.split('by size:')[1].split()}
out = subprocess.run([sys.executable, os.path.join(here, 'closed_sets_ext.py'), '3', str(CMAX)], check=True, capture_output=True, text=True).stdout
m = re.findall(r'\{[^{}]*\}', out)
cb = {int(k): int(v) for k, v in re.findall(r'(\d+): (\d+)', m[-1])} if m else {}
same = {n: c for n, c in ca.items() if c} == {n: c for n, c in cb.items() if c}; ok = ok and same
print(f"Z^3, closed polymers through a plaquette, n <= {CMAX}: {ca} (closed_sets.c) == {cb} (closed_sets_ext.py): {same}", flush=True)
closed = {'3': {str(n): ca.get(n, 0) for n in range(10, CMAX + 1)}}
if not ok: raise SystemExit('the two programs disagree')
J = {'rooted': rooted, 'oneplaq': oneplaq, 'closed': closed,
     'note': "Lemma mg:lem:smallcounts. rooted: {d: {'4' or '6': {n: S_n}}} in Z^(d-1). oneplaq: {d: {E0, terms: [count, 2E, |p'|, m], empty: classes with no junction}}"}
json.dump(J, open(os.path.join(here, 'orders', 'sharp_counts.json'), 'w'), indent=1)
print('all counts made by two programs agree; written to orders/sharp_counts.json')
