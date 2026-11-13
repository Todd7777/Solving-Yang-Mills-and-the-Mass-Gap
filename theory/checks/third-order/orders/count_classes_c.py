# Section 16 of Part IV: junction classes counted by junctions.c (the search of order_counts.py in C) and stored in a cache in the format of
# order_counts_par.py ("[k1, k2, ...]|nmax" -> {n: count}).
# Usage: python3 orders/count_classes_c.py D cache.json JUNCTIONS_BINARY K1:NMAX K2:NMAX ...     (K as k1,k2,...; e.g. 4,4:7 4,6:6)
"""Each class (K, n_max) is counted once (classes already in the cache are skipped), smallest first, and the cache is written after each one, so an
interrupted run resumes. The binary is built from junctions.c (cc -O2 -o junctions junctions.c)."""
import sys, json, subprocess, re, time

def main():
    D, cache_path, binary = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    todo = []
    for arg in sys.argv[4:]:
        ks, nmax = arg.split(':')
        K = sorted(int(x) for x in ks.split(','))
        todo.append((K, int(nmax)))
    try: cache = json.load(open(cache_path))
    except Exception: cache = {}
    todo.sort(key=lambda t: (max(t[0]), t[1], len(t[0])))
    for K, nmax in todo:
        key = f"{K}|{nmax}"
        if key in cache: print(f"{key}: cached {cache[key]}"); continue
        t0 = time.time()
        out = subprocess.run([binary, str(D), str(nmax), ','.join(map(str, K))], capture_output=True, text=True, check=True).stdout
        m = re.search(r"junctions by n:(.*)", out)
        cnt = {} if m is None else {p.split(':')[0]: int(p.split(':')[1]) for p in m.group(1).split()}
        cache[key] = cnt
        json.dump(cache, open(cache_path, 'w'), indent=0)
        print(f"{key}: {cnt}  ({time.time() - t0:.0f} s)"); sys.stdout.flush()

if __name__ == "__main__":
    main()
