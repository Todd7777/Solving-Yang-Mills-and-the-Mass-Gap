# Referee (Section 16, seventh pass): the numbers of Lemma mg:lem:frontsdecision ("Fronts counted by decisions").
# Usage: python3 decision_lemma_check.py
"""(1) The constants quoted: C(s) = (1+s)^{2D} - 1 - 2Ds and y(s)/s, y = max((1+s)^{2D-1} - 1, (1+s)^{D-1}), for D = 3, s = 0.194 and D = 2,
s = 0.325, and the s minimising the growth y/s.  (2) The bound against the exact N^v_k (referee's fronts_redel.c; Z^3 to k = 14, Z^2 to 16 from the
sixth referee's enumeration).  (3) The potential inequality of the proof, checked by brute force at a single decision: for every degree deg_S(w) in
{1, 2, ..., 2D-1} and every number u of undecided edges at w (u <= 2D - deg), the sum over admissible T of s^{|T|} is at most the factor
y^{min(deg, 2)} that the decision removes.  (4) The rate at which the barred tree sums converge with the bound, 7.355 sqrt(abar) e^{xi_F}, at the
abar = 2^-11 chosen at beta_0(4), against 25.8 sqrt(abar) for (d-1) binom((4d-6)k, k-1)."""
from mpmath import mp, mpf, binomial, sqrt, log, findroot
mp.dps = 30

def C(D, s): return (1 + s)**(2*D) - 1 - 2*D*s
def y(D, s): return max((1 + s)**(2*D - 1) - 1, (1 + s)**(D - 1))

for D, s, Cq, gq in ((3, mpf('0.194'), mpf('0.734'), mpf('7.355')), (2, mpf('0.325'), mpf('0.783'), mpf('4.081'))):
    print(f"D={D}, s={s}: C = {float(C(D, s)):.6f} (quoted {Cq}: {'ok' if C(D, s) <= Cq else 'TOO SMALL'}), y = {float(y(D, s)):.6f}, "
          f"y/s = {float(y(D, s)/s):.6f} (quoted {gq}: {'ok' if y(D, s)/s <= gq else 'TOO SMALL'})")
    # the s at which the two arms of y are equal minimises y/s on the right arm's side; scan
    best = min((y(D, mpf(i)/10000)/(mpf(i)/10000), mpf(i)/10000) for i in range(100, 10000))
    print(f"   min over s of y/s: {float(best[0]):.6f} at s = {float(best[1]):.4f}")

EXACT3 = {4: 12, 6: 132, 7: 108, 8: 2034, 9: 4344, 10: 43650, 11: 147888, 12: 1121083, 13: 4773186, 14: 31228116}
EXACT2 = {4: 4, 6: 12, 7: 12, 8: 70, 9: 160, 10: 634, 11: 1728, 12: 6369, 13: 18704, 14: 65228, 15: 199224, 16: 671188}
for D, ex, s in ((3, EXACT3, mpf('0.194')), (2, EXACT2, mpf('0.325'))):
    ok = all(v <= C(D, s)*(y(D, s)/s)**k for k, v in ex.items())
    tot = sum(v*(s/y(D, s))**k for k, v in ex.items())
    ratios = ', '.join(f"{k}:{float(v/(C(D, s)*(y(D, s)/s)**k)):.1e}" for k, v in sorted(ex.items()))
    print(f"D={D}: exact N^v_k <= bound for all listed k: {ok}; sum_k N^v_k (s/y)^k over the listed k = {float(tot):.4f} <= C = {float(C(D, s)):.4f};"
          f" exact/bound: {ratios}")

# (3) one decision: deg = deg_S(w) >= 1, u undecided edges (u <= 2D - deg); T ranges over subsets of the u edges with deg + |T| >= 2
for D, s in ((3, mpf('0.194')), (2, mpf('0.325')), (3, mpf('0.05')), (3, mpf('0.5')), (4, mpf('0.15'))):
    Y = y(D, s); worst = 0
    for deg in range(1, 2*D):
        for u in range(0, 2*D - deg + 1):
            tot = sum(binomial(u, t)*s**t for t in range(0, u + 1) if deg + t >= 2)
            worst = max(worst, tot/Y**min(deg, 2))
    print(f"D={D}, s={s}: max over (deg, u) of sum_T s^|T| / y^min(deg,2) = {float(worst):.6f} (<= 1 up to rounding: {worst <= 1 + mpf(10)**-25})")

# (4) convergence rates of the barred sums sum_k N^v_k abar^{k/2} e^{xi_F k}
abar = mpf(2)**-11
print(f"abar = 2^-11: sqrt(abar) = {float(sqrt(abar)):.5f};  7.355 sqrt(abar) = {float(7.355*sqrt(abar)):.4f};  "
      f"25.81 sqrt(abar) = {float(mpf('25.81')*sqrt(abar)):.4f};  binom(10k,k-1)^(1/k) at k = 60: {float(binomial(600, 59)**(mpf(1)/60)):.3f}")
