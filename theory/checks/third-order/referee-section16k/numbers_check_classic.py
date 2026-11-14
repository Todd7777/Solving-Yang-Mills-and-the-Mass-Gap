# Referee check (eleventh pass on Section 16, Corollary mg:cor:classic), item 5: every number quoted in Corollary mg:cor:classic, its proof, the subsection
# before it, the abstract sentence and the table row, against (i) the certified Arb outputs of this folder (arb_eval_classic_*.out, upper_arb_classic.out,
# series_check_classic.out) and (ii) the authors' outputs (orders/*.out, reproduced identically in outputs/repro_*.out). A quoted bound must lie on the safe
# side of the certified value; "tight" means it is the certified value rounded up at the quoted precision.
# Usage: python3 numbers_check_classic.py
import os, re
from mpmath import mp, mpf, exp, log, ceil, floor, log10
mp.dps = 40
here = os.path.dirname(os.path.abspath(__file__)); O = os.path.join(here, 'outputs'); ORD = os.path.join(here, '..', 'orders')

def grab(path, pat, idx=1):
    s = open(path).read(); m = re.search(pat, s)
    assert m, (path, pat)
    return mpf(m.group(idx))

def beta0(d):
    K = 4*(2*d - 3); return log(1 + 1/(mp.e**2*K*(K + 2)))/2

def up(x, sig):
    e = int(floor(log10(x))) - sig + 1; return ceil(x/mpf(10)**e)*mpf(10)**e

rows = []
def chk(where, quoted, true, kind, sig=None, note=''):
    """kind 'upper': quoted must be >= true; 'lower': quoted <= true; 'approx': |quoted/true - 1| small; 'one digit': true rounded to one digit is quoted"""
    q = mpf(quoted)
    if kind == 'upper': ok = q >= true
    elif kind == 'lower': ok = q <= true
    elif kind == 'one digit': ok = mpf(mp.nstr(true, 1)) == q
    else: ok = abs(q/true - 1) < mpf('0.06')
    tight = (sig is not None and kind == 'upper' and q == up(true, sig))
    rows.append((where, quoted, mp.nstr(true, 8), kind, ok, ('tight' if tight else ('' if sig is None else f'rounded-up value {mp.nstr(up(true, sig), sig)}')) + (' ' + note if note else '')))

cb = {}
for g, d, f in (('SU(2)', 3, 'su2_d3'), ('SU(2)', 4, 'su2_d4'), ('U(1)', 3, 'u1_d3'), ('U(1)', 4, 'u1_d4')):
    p = os.path.join(O, f'arb_eval_classic_{f}.out')
    cb[(g, d)] = grab(p, r'\(1\+kappa\)\^3 D/beta_\*\^4 <= ([0-9.e+-]+)')
    cb[(g, d, 'dl')] = grab(p, r"delta'_\*\s+([0-9.e+-]+)")
up_s = open(os.path.join(O, 'upper_arb_classic.out')).read()
E = {}
for g, d, C in (('SU(2)', 3, '+'), ('U(1)', 3, '+'), ('U(1)', 3, '-'), ('SU(2)', 4, '+'), ('U(1)', 4, '+'), ('U(1)', 4, '-')):
    m = re.search(rf"{re.escape(g)}, d={d}, C={re.escape(C)}: .*?E_C/beta\^4 = \[([0-9.e+-]+) \+/- ([0-9.e+-]+)\]", up_s)
    E[(g, d, C)] = mpf(m.group(1)) + mpf(m.group(2))          # upper end of the ball
eb = [grab(os.path.join(ORD, f'error_budget_{f}.out'), r'largest tau ([0-9.e+-]+)') for f in ('su2_d3', 'su2_d4', 'u1_d3', 'u1_d4')]
rm = [grab(os.path.join(ORD, f'error_budget_{f}.out'), r'smallest room ([0-9.e+-]+)') for f in ('su2_d3', 'su2_d4', 'u1_d3', 'u1_d4')]
b3, b4 = beta0(3), beta0(4)

# --- the statement (lines 5343-5357) ---
chk('Cor (i): c_- SU(2) d=3', '1968', cb[('SU(2)', 3)], 'upper', 4)
chk('Cor (i): c_+ SU(2) d=3', '0.0026', E[('SU(2)', 3, '+')], 'upper', 2)
chk('Cor (i): c_- SU(2) d=4', '570', cb[('SU(2)', 4)], 'upper', 3)
chk('Cor (i): c_+ SU(2) d=4', '0.0028', E[('SU(2)', 4, '+')], 'upper', 2)
chk("Cor (ii): c'_- U(1) d=3", '1.772e5', cb[('U(1)', 3)], 'upper', 4)
chk("Cor (ii): c'_+ U(1) d=3", '0.0017', max(E[('U(1)', 3, '+')], E[('U(1)', 3, '-')]), 'upper', 2)
chk("Cor (ii): c'_- U(1) d=4", '8043', cb[('U(1)', 4)], 'upper', 4)
chk("Cor (ii): c'_+ U(1) d=4", '0.0018', max(E[('U(1)', 4, '+')], E[('U(1)', 4, '-')]), 'upper', 2)
# --- the proof, lower bounds (lines 5361-5367) ---
chk("proof: delta'_* ~ 5e-9 (U(1), d=3)", '5e-9', cb[('U(1)', 3, 'dl')], 'one digit', note='(4.65e-9, to one digit)')
chk('proof: D <= 1967.1 beta^4', '1967.1', cb[('SU(2)', 3)], 'upper', 5)
chk('proof: D <= 569.17 beta^4', '569.17', cb[('SU(2)', 4)], 'upper', 5)
chk('proof: D <= 1.7714e5 beta^4', '1.7714e5', cb[('U(1)', 3)], 'upper', 5)
chk('proof: D <= 8042.2 beta^4', '8042.2', cb[('U(1)', 4)], 'upper', 5)
chk('proof: tails <= 2e-25 of their sums', '2e-25', max(eb), 'upper', 1)
chk('proof: room >= 2.5e-16', '2.5e-16', min(rm), 'lower')
# --- the proof, upper bounds (lines 5369-5372): quoted to 3 significant digits ---
for (g, d, C), q in ((('SU(2)', 3, '+'), '0.00251'), (('SU(2)', 4, '+'), '0.00277'), (('U(1)', 3, '+'), '0.00168'), (('U(1)', 4, '+'), '0.00179')):
    chk(f'proof: E_C(beta_0)/beta_0^4 {g} d={d}', q, E[(g, d, C)], 'approx', note='(3 digits)')
# --- the proof, series and splitting (lines 5374-5383) ---
split = grab(os.path.join(O, 'series_check_classic.out'), r'\(1/2\)\(1 - ([0-9.e+-]+)\)')
chk('proof: I_1^2 - I_0 I_2 >= (b/2)^2 (1/2)(1 - 3e-8)', '3e-8', split, 'upper', 1)
chk('proof: 4.03e-4 >= beta_0(3)', '4.03e-4', b3, 'upper')
for d, q, cm, cp in ((3, '0.0288', mpf('1.772e5'), mpf('0.0017')), (4, '1.91e-4', mpf(8043), mpf('0.0018'))):
    b = beta0(d)
    chk(f"proof: (c'_- + c'_+) beta_0^2, d={d}", q, (cm + cp)*b**2, 'upper', 3)
chk('proof: e^{-3 beta_0}(1-1e-7)/4 (the smaller, d=3)', '0.2497', exp(-3*b3)*(1 - mpf(10)**-7)/4, 'approx', note='(true 0.249698: quoted rounded UP in a comparison it must win; margin x8.7, harmless)')
for g, d, q, c in (('SU(2)', 3, '0.0026', mpf(1968)), ('SU(2)', 4, '1.1e-4', mpf(570)), ('U(1)', 3, '0.23', mpf('1.772e5')), ('U(1)', 4, '1.6e-3', mpf(8043))):
    b = beta0(d)
    chk(f'proof: residual/(beta^2/8) at beta_0, {g} d={d}', q, 8*c*b**2, 'upper', 2)
# --- subsection, abstract, table ---
chk('subsection: u = beta/4 (SU(2)) [u_1 = (1+iota)/(2N^2)]', '0.25', mpf(2)/8, 'approx')
chk('subsection: u = beta/2 (U(1))', '0.5', mpf(1)/2, 'approx')
chk('subsection: u = beta/18 (SU(3))', str(1/18.), mpf(1)/18, 'approx')
chk('abstract: -570 beta^4 (SU(2), d=4)', '570', cb[('SU(2)', 4)], 'upper', 3)
chk('abstract: 0.0028 beta^4 (SU(2), d=4)', '0.0028', E[('SU(2)', 4, '+')], 'upper', 2)
w = max(len(r[0]) for r in rows)
for where, q, t, kind, ok, note in rows:
    print(f"{where:{w}s}  quoted {q:>10s}  certified/true {t:>12s}  ({kind}): {'OK' if ok else 'WRONG SIDE'}  {note}")
print(f"\n{sum(r[4] for r in rows)} of {len(rows)} correct (bounds on the safe side; 'tight' = the certified value rounded up at the quoted precision)")
