# Section 16 of Part IV: junctions.c (the search of order_counts.py in C) against the counts of the Python searches.
# Usage: python3 junctions_validate.py JUNCTIONS_BINARY        (build: cc -O2 -o junctions junctions.c; about 3 minutes)
"""Every class below was counted by order_counts.py (count_junctions, Python) or junction_count.py; the C search must reproduce each count, for every
n <= n_max (an absent n means 0). Sources: orders/d3_m8_classes.txt and orders/assemble_d3_s9.py (d=3), orders/d4_m8_classes.json and Remark
mg:rem:practical(ii), (ix) (d=4); the d=4 classes of order 7 are those of orders/d4_m7.json (order_counts_par.py 3 7)."""
import sys, subprocess, re

CASES = [  # (D, n_max, sizes, {n: count}) ; D = d - 1
    # d=3, order 8 (orders/d3_m8_classes.txt) and order 9 (orders/assemble_d3_s9.py)
    (2, 6, (4, 4), {3: 48}), (2, 7, (4, 4), {3: 48}),
    (2, 6, (4, 6), {2: 24, 3: 24, 4: 276}), (2, 6, (4, 7), {2: 24, 3: 24, 4: 276}),
    (2, 5, (4, 8), {1: 4, 2: 48, 3: 252, 4: 604, 5: 1152}), (2, 5, (4, 9), {1: 8, 2: 16, 3: 424, 4: 704, 5: 2288}),
    (2, 4, (4, 10), {1: 20, 2: 100, 3: 732, 4: 2432}), (2, 4, (4, 11), {1: 16, 2: 212, 3: 856, 4: 4556}),
    (2, 3, (4, 12), {1: 12, 2: 380, 3: 1872}), (2, 2, (4, 13), {2: 412}),
    (2, 4, (6, 6), {3: 64, 4: 64}), (2, 4, (6, 7), {3: 128, 4: 128}), (2, 2, (6, 8), {2: 48}), (2, 2, (6, 9), {2: 60}), (2, 2, (7, 7), {}),
    (2, 5, (4, 4, 4), {4: 292, 5: 4}), (2, 4, (4, 4, 6), {3: 168, 4: 168}), (2, 4, (4, 4, 7), {3: 168, 4: 168}),
    (2, 2, (4, 4, 8), {}), (2, 2, (4, 4, 9), {}), (2, 2, (4, 6, 6), {2: 20}), (2, 2, (4, 6, 7), {2: 40}), (2, 2, (4, 4, 4, 4), {}),
    # d=3: Z_(n,4) (junction_count.py: Z_1 = 1, Z_2 = 8, none with 3..7)
    (2, 7, (4,), {1: 1, 2: 8}),
    # d=4: Z_(n,4) to 9 (junction_count.py 3 9)
    (3, 9, (4,), {1: 1, 2: 32, 4: 12, 5: 26, 6: 302, 7: 288, 8: 1988, 9: 406}),
    # d=4, the classes of order 8 counted in Python (orders/d4_m8_classes.json) and of order 7
    (3, 6, (4, 4), {3: 876, 4: 42, 5: 912, 6: 1866}), (3, 5, (4, 6), {2: 320, 3: 680, 4: 18276, 5: 11744}),
    (3, 5, (4, 7), {2: 312, 3: 344, 4: 17028, 5: 3440}), (3, 4, (4, 8), {1: 36, 2: 992, 3: 9220, 4: 56960}),
    (3, 4, (4, 9), {1: 72, 2: 256, 3: 15576, 4: 50920}), (3, 3, (4, 10), {1: 276, 2: 2260, 3: 45700}),
    (3, 4, (4, 4, 4), {2: 4, 3: 16, 4: 24314}), (3, 3, (6, 6), {2: 64, 3: 3028}), (3, 3, (6, 7), {3: 5712}),
    (3, 3, (4, 4, 6), {3: 9768}), (3, 3, (4, 4, 7), {3: 9600}), (3, 2, (4, 11), {1: 240, 2: 4636}),
]

def run(binary, D, nmax, sizes):
    out = subprocess.run([binary, str(D), str(nmax), ','.join(map(str, sizes))], capture_output=True, text=True, check=True).stdout
    m = re.search(r"junctions by n:(.*)", out)
    return {} if m is None else {int(p.split(':')[0]): int(p.split(':')[1]) for p in m.group(1).split()}

if __name__ == "__main__":
    binary = sys.argv[1]; bad = 0
    for D, nmax, sizes, want in CASES:
        got = run(binary, D, nmax, sizes)
        ok = got == want
        bad += not ok
        print(f"D={D} n<={nmax} sizes={sizes}: {got} {'ok' if ok else 'DIFFERS from ' + str(want)}"); sys.stdout.flush()
    print(f"{len(CASES)} classes, {bad} differences")
