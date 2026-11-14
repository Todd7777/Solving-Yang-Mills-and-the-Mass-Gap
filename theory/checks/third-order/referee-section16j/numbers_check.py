# Referee check (tenth pass on Section 16, Corollary mg:cor:sectorscap), item 6: every number quoted in the corollary, its proof, the paragraph after it and
# the subsection before it, against (a) the authors' stored outputs in ../orders/ and (b) the referee's ball-arithmetic outputs in outputs/.
# "bound" = the quoted number is used as an upper (lower) bound and must be >= (<=) the certified value; "value" = an approximate value (rounding checked).
# Usage: python3 numbers_check.py       (after arb_eval_cor.py, upper_arb_sectors.py and split_check_arb.py have written their outputs)
import os, re
here = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(here, 'outputs'); A = os.path.join(here, '..', 'orders')
def read(p): return open(p).read()
def num(pat, text, grp=1): m = re.search(pat, text); assert m, pat; return float(m.group(grp))

rows = []
def bound_up(name, quoted, cert):      # quoted must be >= the certified upper bound
    rows.append((name, quoted, cert, 'upper bound', quoted >= cert))
def bound_lo(name, quoted, cert):      # quoted must be <= the certified lower bound
    rows.append((name, quoted, cert, 'lower bound', quoted <= cert))
def value(name, quoted, cert, digits=None):  # quoted is cert rounded to as many significant digits as the quoted literal has
    from math import log10, floor
    from decimal import Decimal
    if digits is None or True:
        digits = len(Decimal(repr(quoted)).normalize().as_tuple().digits)
    r = round(cert, -int(floor(log10(abs(cert)))) + digits - 1)
    rows.append((name, quoted, cert, f'value ({digits} digits)', abs(r - quoted) < 1e-12*abs(cert) + 1e-15))

# lower constants: (1+kappa)^3 D(beta_*)/beta_*^4 in Arb
for (G, d, s, cm, dq) in ((3, 3, 's9z', 0.600, 0.59934), (4, 3, 's9z', 0.134, 0.13324), (5, 3, 's9z', 0.073, 0.072888), (4, 4, 'x92', 1.083, 1.0822),
                          (3, 4, 'x92', 3.937, None), (5, 4, 'x92', 0.520, None)):
    t = read(os.path.join(O, f'arb_eval_cor_su{G}_d{d}_{s}.out'))
    c = num(r'\(1\+kappa\)\^3 D/beta_\*\^4 <= ([0-9.]+)', t)
    assert 'all certified: True' in t
    bound_up(f"c'_- SU({G}) d={d} (table)", cm, c)
    if dq is not None: bound_up(f"proof: D <= {dq} beta^4, SU({G}) d={d}", dq, c)
    # the authors' scaling output, same number
    fn = {('s9z', 3): f'capz_su{G}_d3_scaling.out', ('x92', 4): (f'cap_su{G}_d4_scaling.out')}[(s, d)]
    if os.path.exists(os.path.join(A, fn)):
        ta = read(os.path.join(A, fn)); ca = num(r'\(1\+k\)\^3 D_\*/beta_\*\^4 = ([0-9.]+)', ta)
        value(f"  authors' scaling output SU({G}) d={d} (printed)", ca, c, len(re.sub(r'^0\.0*|\.', '', str(ca))))
# the majorant M_S at beta_*
for (G, s, q) in ((3, 's9', 2.354), (3, 's9z', 0.399), (5, 's9', 0.280), (5, 's9z', 0.007)):
    t = read(os.path.join(O, f'arb_eval_cor_su{G}_d3_{s}.out'))
    m = float(re.search(r'one-plaquette, multi ([^\n]+)', t).group(1).split(',')[-1])
    value(f"proof: majorant M_S SU({G}) d=3 set {s}", q, m, 3 if q > 0.01 else 1)
# upper constants and the quoted E_C(beta_0)/beta_0^4
t = read(os.path.join(O, 'upper_arb_sectors.out'))
quotedE = {(3, 3): 0.04153, (3, 4): 0.04697, (4, 3): 0.13057, (4, 4): 0.14783, (5, 3): 0.31834, (5, 4): 0.36051}
cplus = {(3, 3): 0.042, (3, 4): 0.047, (4, 3): 0.131, (4, 4): 0.148, (5, 3): 0.319, (5, 4): 0.361}
for (G, d) in quotedE:
    vals = [float(re.search(rf'SU\({G}\), d={d}, C={re.escape(C)}:.*?E_C/beta\^4 = \[([0-9.]+) ', t).group(1)) for C in '+-']
    bound_up(f"c'_+ SU({G}) d={d} (table), both sectors", cplus[(G, d)], max(vals))
    value(f"proof: E_C(beta_0)/beta_0^4 SU({G}) d={d}", quotedE[(G, d)], max(vals), 4)
# splitting and the comparisons
t = read(os.path.join(O, 'split_check_arb.out'))
r = re.findall(r'quoted: e\^\(-3b\)\(1-1e-8\)/3000 = \[([0-9.]+) .*?\(c\'_-\+c\'_\+\) beta = \[([0-9.]+) .*?factor \[([0-9.]+) ', t)
(l3, r3, f3), (l4, r4, f4) = [(float(a), float(b), float(c)) for a, b, c in r]
value("proof: e^{-3b}(1-1e-8)/3000 at beta_0(3)", 3.329e-4, l3, 4); value("proof: (c'_-+c'_+) beta_0(3), SU(5)", 1.58e-4, r3, 3)
value("proof: e^{-3b}(1-1e-8)/3000 at beta_0(4)", 3.332e-4, l4, 4); value("proof: (c'_-+c'_+) beta_0(4), SU(5)", 1.36e-4, r4, 3)
bound_lo("after the proof: margin factor d=3 (as a lower bound)", 2.1, f3); value("after the proof: margin factor d=3", 2.1, f3, 2)
bound_lo("after the proof: margin factor d=4 (as a lower bound)", 2.4, f4); value("after the proof: margin factor d=4", 2.4, f4, 2)
r = re.findall(r'against 1/\(648 beta_0\) = \[([0-9.]+) .*?beta\^3/1296 >= \[([0-9.]+) ', t)
(a3, b3), (a4, b4) = [(float(x), float(y)) for x, y in r]
value("proof: 1/(648 beta_0(3))", 3.83, a3, 3); value("proof: 1/(648 beta_0(4))", 10.04, a4, 4)
bound_lo("(i): beta^3/1296 >= 1.91 beta^4 (d=3)", 1.91, b3); bound_lo("(i): beta^3/1296 >= 5.01 beta^4 (d=4)", 5.01, b4)
m = re.findall(r'margin = \[([0-9.]+) .*?q_N \(1 - \[([0-9.e-]+) ', t)
for (G, (mg, rel)) in zip((3, 4, 5), m): bound_up(f"proof: Q >= q_N t^(N-2) (1 - 1e-8), SU({G}): relative loss", 1e-8, float(rel))
# error budget of the new cases
t = read(os.path.join(A, 'error_budget_sectors.out'))
taus, rooms = [], []
for line in t.splitlines():
    mm = re.match(r'SU\((\d)\), d=(\d)( \(d3_s9z.json\))?: .*largest relative tail tau = ([0-9.e-]+) .*?room with the tails included: (.*?);  all', line)
    if not mm: continue
    G, d = int(mm.group(1)), int(mm.group(2))
    if (G, d) in ((3, 4), (5, 4)): continue                  # not new cases
    taus.append(float(mm.group(4))); rooms.append(min(float(x.split()[-1]) for x in mm.group(5).split(', ')))
bound_up("proof: tails at most 3e-27 of their sums (new cases)", 3e-27, max(taus)); bound_lo("proof: room at least 2.8e-16 (new cases)", 2.8e-16, min(rooms))
w = max(len(r[0]) for r in rows)
bad = 0
for name, q, c, kind, ok in rows:
    bad += not ok
    print(f"{name:{w}s}  quoted {q:<10g} computed {c:<14.10g} {kind:18s} {'ok' if ok else 'DIFFERS'}")
print(f"{len(rows)} numbers, {bad} flagged")
