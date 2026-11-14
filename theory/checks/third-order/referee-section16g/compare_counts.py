# Referee (Section 16, seventh pass): the junction counts of the partition search (junc_partition_big.c, on the referee's own front dump) against the
# authors' class caches (orders/d4_m8_classes.json, d4_m9_classes.json, d4_m10_classes.json) and Z_(n,4) (practical_beta.Z4_D4).
# Usage: python3 compare_counts.py OUT1.txt [OUT2.txt ...]
import json, os, re, sys
here = os.path.dirname(os.path.abspath(__file__)); orders = os.path.join(here, '..', 'orders')
auth = {}
for f in ('d4_m8_classes.json', 'd4_m9_classes.json', 'd4_m10_classes.json'):
    for key, cnt in json.load(open(os.path.join(orders, f))).items():
        K = tuple(json.loads(key.split('|')[0])); nm = int(key.split('|')[1])
        if nm >= auth.get(K, (-1, None))[0]: auth[K] = (nm, {int(n): v for n, v in cnt.items()})
auth[(4,)] = (12, {1: 1, 2: 32, 4: 12, 5: 26, 6: 302, 7: 288, 8: 1988, 9: 406, 10: 3193, 11: 5911, 12: 26956})   # Z_(n,4): Remark 16.5(ii), (ix), (x)
rows = []; bad = 0
for path in sys.argv[1:]:
    txt = open(path).read()
    for m in re.finditer(r"NMAX=(\d+) sizes=([\d,]+) prune=\d.*?\n(?:  junctions by n:(.*?)\n|hash set full)", txt, re.S):
        nmax, K = int(m.group(1)), tuple(sorted(int(x) for x in m.group(2).split(',')))
        if m.group(3) is None: rows.append((K, nmax, 'hash set full', '')); continue
        got = {int(p.split(':')[0]): int(p.split(':')[1]) for p in m.group(3).split()}
        if K not in auth: rows.append((K, nmax, got, 'not in the authors\' caches')); continue
        anm, ac = auth[K]
        cmp = {n: (got.get(n, 0), ac.get(n, 0)) for n in range(1, min(nmax, anm) + 1)}
        diff = {n: v for n, v in cmp.items() if v[0] != v[1]}
        bad += bool(diff)
        rows.append((K, nmax, {n: v for n, v in got.items() if v}, 'equal' if not diff else f'DIFFERS {diff}'))
for K, nmax, got, verdict in rows:
    print(f"K = {list(K)}, n <= {nmax}: referee {got}  -> {verdict}")
print(f"{len(rows)} runs, {bad} differences")
