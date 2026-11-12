# Written by the independent referee of Section 15 of Part IV (2026-11-12) and added to the repository with credit; unchanged apart from
# this header. A compass, not part of the proofs.
# Referee of Section 15 (2026-11-12). A compass, not part of the proofs.
# Corollary 15.10(i) for U(1) (rho = e^{i theta}): Upsilon_- = 1 - I_2/I_0 (alpha = 0, gamma = I_2/I_0), Upsilon_+ = 1 + I_2/I_0 - 2 (I_1/I_0)^2,
# u = I_1/I_0. Checks: log Upsilon_- - log Upsilon_+ = beta^2/4 + O(beta^3); -log Upsilon_- = u^2/2 + (7/24) u^4 + O(u^6), so with h_- = 10d-22
# the u^4 coefficient of m + 4 log u is 7/24 - 18 = -425/24 (d=4) and 7/24 - 8 = -185/24 (d=3); and R(beta) -> -h_- u_1^4 = -(10d-22)/16.
import mpmath as mp
mp.mp.dps = 60
for beta in (mp.mpf('1e-2'), mp.mpf('1e-3'), mp.mpf('1e-4')):
    I0, I1, I2 = (mp.besseli(k, beta) for k in (0, 1, 2))
    u = I1 / I0
    Um = 1 - I2 / I0
    Up = 1 + I2 / I0 - 2 * u**2
    L = -mp.log(Um)
    print(f'beta={mp.nstr(beta,3)}: (log U- - log U+)/beta^2 = {mp.nstr((mp.log(Um)-mp.log(Up))/beta**2, 12)}  (1/4)'
          f'   (-log U- - u^2/2)/u^4 = {mp.nstr((L - u**2/2)/u**4, 12)}  (7/24 = {mp.nstr(mp.mpf(7)/24, 12)})')
print('u^4 coefficients of m + 4 log u (C = -):', 'd=4:', mp.mpf(7)/24 - 18, '= -425/24 =', mp.mpf(-425)/24, '; d=3:', mp.mpf(7)/24 - 8, '= -185/24 =', mp.mpf(-185)/24)
print('R(beta) limit -(10d-22)/16 with u_1 = 1/2: d=4:', -(10*4-22)/16, ' d=3:', -(10*3-22)/16)
