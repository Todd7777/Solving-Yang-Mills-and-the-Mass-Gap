# Section 16 of Part IV: order_counts.py with the classes run in parallel and cached.  Usage: python3 order_counts_par.py D m JOBS cache.json out.json
"""The same counts and terms as order_counts.py; each class (K, n_max) is counted once by count_junctions and stored in the cache file, so an
interrupted run resumes.  The classes are taken longest-first."""
import sys, json, multiprocessing as mp
from collections import Counter
from fractions import Fraction
import order_counts as OC
def work(args):
    D, K, nmax = args
    c = OC.count_junctions(D, nmax, K)
    return (list(K), nmax, {str(n): v for n, v in c.items()})
if __name__ == "__main__":
    D, m, jobs, cache_path, out_path = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
    try: cache = json.load(open(cache_path))
    except Exception: cache = {}
    C = OC.classes(m)
    todo = [(D, K, nmax) for K, nmax in C.items() if f"{list(K)}|{nmax}" not in cache]
    todo.sort(key=lambda t: (-max(t[1]), -t[2]))
    print(f"{len(C)} classes, {len(todo)} to count"); sys.stdout.flush()
    with mp.Pool(jobs) as pool:
        for K, nmax, cnt in pool.imap_unordered(work, todo):
            cache[f"{K}|{nmax}"] = cnt
            json.dump(cache, open(cache_path, 'w'))
            print(f"K={tuple(K)} (n<={nmax}): {cnt}"); sys.stdout.flush()
    lim = Fraction(2 * m + 1, 2); agg = Counter()
    for K, nmax in C.items():
        cnt = {int(n): v for n, v in cache[f"{list(K)}|{nmax}"].items()}
        for n in range(2, nmax + 1):
            N = cnt.get(n, 0)
            if not N: continue
            for k, mult in Counter(K).items():
                Ek = OC.E_of(K, k)
                if n + Ek <= lim and (n > 2 or OC.E_of(K, K[0]) >= 3): agg[(n, float(Ek), k, len(K) - 1)] += mult * N
    terms = [(c, n, Ek, k, r1) for (n, Ek, k, r1), c in sorted(agg.items())]
    corners = [(n, m + 1 - n) for n in range(2, m) if (m + 1 - n) >= (3 if n == 2 else 2)]
    json.dump({'terms': terms, 'corners': corners}, open(out_path, 'w'))
    print("terms:", len(terms), "corners:", corners)
