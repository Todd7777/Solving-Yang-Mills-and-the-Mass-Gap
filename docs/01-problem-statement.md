# 1. The problem, stated precisely

This document fixes *exactly* what has to be proved. Every later document in this repository
refers back to the definitions here. Nothing in this file is new mathematics; its purpose is to
remove ambiguity so that partial results can be measured against the real target.

## 1.1 The Clay Mathematics Institute statement

The official problem description (A. Jaffe and E. Witten, *Quantum Yang–Mills Theory*, 2000)
asks:

> Prove that for any compact simple gauge group $G$, a non-trivial quantum Yang–Mills theory
> exists on $\mathbb{R}^4$ and has a mass gap $\Delta > 0$. Existence includes establishing
> axiomatic properties at least as strong as those cited in [Streater–Wightman; Osterwalder–Schrader].

Three things are being asked for, and each is a separate mathematical statement:

1. **Existence.** A quantum field theory on four-dimensional Minkowski space (equivalently, by
   Osterwalder–Schrader reconstruction, a Euclidean theory on $\mathbb{R}^4$) satisfying the
   Wightman axioms (equivalently, the Osterwalder–Schrader axioms), whose local observables are the
   gauge-invariant composites of a $G$-connection.
2. **Non-triviality.** The theory is not a (generalized) free field. In practice one wants more:
   the theory should be *the* theory whose short-distance behaviour is governed by the classical
   Yang–Mills action with asymptotic freedom. The problem statement only demands non-triviality.
3. **Mass gap.** The Hamiltonian $H \ge 0$ has a simple eigenvalue $0$ (the vacuum $\Omega$) and
   $\operatorname{spec}(H) \subset \{0\} \cup [\Delta, \infty)$ for some $\Delta > 0$.

## 1.2 Classical Yang–Mills

Let $G$ be a compact simple Lie group with Lie algebra $\mathfrak{g}$, equipped with a faithful
unitary representation (for $G = SU(N)$, the defining representation) and the invariant inner product
$\langle X, Y\rangle = -2\operatorname{Tr}(XY)$ normalised so that $\operatorname{Tr}(T^aT^b)=\tfrac12\delta^{ab}$
for an orthonormal basis $T^a$ of $\mathfrak{g}$ (anti-Hermitian generators absorbed into $i$).
A connection is a $\mathfrak{g}$-valued one-form $A = A_\mu^a T^a\,dx^\mu$ on $\mathbb{R}^4$ (for the
trivial bundle); its curvature is
$$F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu + [A_\mu, A_\nu],$$
and the Euclidean action is
$$S[A] = \frac{1}{2g^2}\int_{\mathbb{R}^4} \operatorname{Tr}(F_{\mu\nu}F_{\mu\nu})\,d^4x
      = \frac{1}{4g^2}\int_{\mathbb{R}^4} F^a_{\mu\nu}F^a_{\mu\nu}\,d^4x .$$
Gauge transformations $A \mapsto gAg^{-1} - (dg)g^{-1}$ leave $S$ invariant. The *gauge-invariant
local observables* are functions of $A$ invariant under gauge transformations: traces of products of
curvature components and their covariant derivatives at a point (e.g. $\operatorname{Tr}F_{\mu\nu}F_{\mu\nu}$),
and Wilson loops $W_\gamma = \operatorname{Tr}\,\mathcal{P}\exp\oint_\gamma A$ for closed curves $\gamma$.

The coupling $g$ is dimensionless in four dimensions. Classically there is no scale; the quantum
theory is expected to generate one by *dimensional transmutation*: the physical mass gap should be
$\Delta = C\,\Lambda$, where $\Lambda$ is the dynamically generated scale and $C$ a pure number
depending only on $G$.

## 1.3 The axioms we target

### Euclidean (Osterwalder–Schrader) axioms

A Euclidean theory is specified by Schwinger functions $S_n(x_1,\dots,x_n)$ for a family of
local fields $\{\Phi_\alpha\}$. We use the standard list (Osterwalder–Schrader 1973/1975;
Glimm–Jaffe, *Quantum Physics*, Chapter 6):

* **(OS0) Temperedness / regularity.** $S_n$ are tempered distributions on $\mathbb{R}^{4n}$ (and,
  for the reconstruction theorem in the form needed to produce Wightman fields, satisfy the
  linear-growth condition $|S_n(f_1\otimes\cdots\otimes f_n)| \le C^n\, n!^{\,\sigma}\prod_i\|f_i\|$
  for some Schwartz norm $\|\cdot\|$ and $\sigma \ge 0$).
* **(OS1) Euclidean invariance.** $S_n$ are invariant under the Euclidean group $E(4)$ (rotations,
  translations; for scalar gauge-invariant fields no spin indices need to be rotated).
* **(OS2) Reflection positivity.** For $\theta$ the reflection $x_0\mapsto -x_0$ and for every finite
  family of test functions $f_i$ supported at positive times, $\sum_{i,j} S(\overline{\theta f_i}\otimes f_j) \ge 0$.
* **(OS3) Symmetry.** $S_n$ are symmetric under permutations of their (bosonic) arguments.
* **(OS4) Clustering.** $S_{n+m}(f\otimes g_{\lambda a}) \to S_n(f)S_m(g)$ as $\lambda\to\infty$ for
  translations $a$ in a spatial direction (ergodicity of the vacuum).

The Osterwalder–Schrader reconstruction theorem turns (OS0)–(OS4) into a Hilbert space $\mathcal{H}$, a
unique vacuum $\Omega$, a self-adjoint $H\ge 0$ with $H\Omega=0$, a positive-energy unitary representation
of the Poincaré group, and Wightman fields.

### Mass gap

Within the reconstructed theory, the **mass gap** is the statement
$$\operatorname{spec}(H)\cap(0,\Delta)=\emptyset \quad\text{for some }\Delta>0,$$
with $0$ a simple eigenvalue. By the spectral theorem this is equivalent to **exponential clustering in
Euclidean time**: for all vectors of the form $A\Omega$, $B\Omega$ with $A,B$ in the algebra generated by
the Euclidean fields at positive times,
$$\bigl|\langle \Omega, A\,e^{-tH}B\,\Omega\rangle-\langle\Omega,A\Omega\rangle\langle\Omega,B\Omega\rangle\bigr|
\le C_{A,B}\,e^{-\Delta t}\qquad(t\ge 0).$$
(Theorem 4.3 of the foundations monograph proves this equivalence in the lattice setting, where it is
needed; the argument is identical in the continuum.)

### What "gauge theory" adds

For a gauge theory the fields $\Phi_\alpha$ are the gauge-invariant composites. The Hilbert space obtained
from gauge-invariant Schwinger functions is automatically the *physical* (gauge-invariant) state space:
no Gupta–Bleuler or BRST subsidiary condition is needed. This is a real advantage of the Euclidean
lattice route and the reason every construction in this project works with gauge-invariant observables only.

## 1.4 The lattice reduction

The only known framework in which all of (OS0)–(OS4) can be approached simultaneously for a
non-abelian gauge theory is the Wilson lattice regularisation (definitions in
`theory/lattice-foundations.tex`, §2). For lattice spacing $a>0$, a finite periodic lattice
$\Lambda_{a,L}$ of physical side $L$ and inverse coupling $\beta=2N/g^2$, Wilson's measure
$d\mu_{a,L,\beta}$ is a well-defined probability measure on $G^{\text{links}}$; it is reflection
positive (Theorem 3.2 of the monograph), has a positive self-adjoint transfer matrix and Hamiltonian
$H_{a,L}$ with a unique vacuum, and a strictly positive gap $m_{a,L}(\beta)>0$ (Theorem 4.2).
The gap is **not** known to be uniform in $L$ or $a$. The problem is to prove it is, along the correct
trajectory $\beta=\beta(a)$.

**Theorem (reduction; folklore, proved in the monograph §5).** Suppose that for a family of
renormalised gauge-invariant local lattice observables $\{\Phi^{(a)}_\alpha\}$ and a function
$\beta(a)\to\infty$ as $a\to 0$:

* **(L1) Convergence.** The lattice Schwinger functions
  $S^{(a,L)}_n$ converge, as $a\to0$ and $L\to\infty$ along $\beta=\beta(a)$, to tempered distributions
  $S_n$ on $\mathbb{R}^{4n}$ satisfying (OS0), and the limits are invariant under the full Euclidean group
  (restoration of rotation invariance from the hypercubic symmetry of the lattice).
* **(L2) Uniform exponential clustering.** There are $\Delta>0$ and, for every pair of test functions
  $f,g$ with bounded supports, a constant $C_{f,g}$ independent of $(a,L)$ with
  $$\bigl|S^{(a,L)}(f\otimes g_t)-S^{(a,L)}(f)\,S^{(a,L)}(g)\bigr|\le C_{f,g}\,e^{-\Delta t}$$
  for all $t\ge0$, where $g_t$ is $g$ translated by $t$ in Euclidean time.
* **(L3) Non-triviality.** Some truncated four-point function of the limit does not vanish identically.

Then the limiting Schwinger functions satisfy (OS0)–(OS4), the reconstructed theory is non-trivial, and its
Hamiltonian has mass gap at least $\Delta$.

The proof is short because (OS2) and (OS3) are inherited from the lattice (reflection positivity passes
to limits of positive quadratic forms), (OS4) follows from (L2), and the gap follows from (L2) by the
spectral argument above. Everything hard is in (L1) and (L2).

## 1.5 Why this is hard: two regimes and a crossover

The hypotheses (L1) and (L2) live at opposite ends of the renormalisation group flow.

* **Ultraviolet (L1).** As $a\to0$, asymptotic freedom forces $\beta(a)\to\infty$: the lattice theory is
  weakly coupled at the cutoff scale. The control of the limit $a\to0$ in finite volume is the content of
  Balaban's renormalisation-group program (1982–1989) and Magnen–Rivasseau–Sénéor (1993), which prove
  *ultraviolet stability* (uniform bounds) but not convergence of the correlation functions.
* **Infrared (L2).** At distances of order $1/\Lambda$ the effective coupling is of order one. No rigorous
  technique controls a non-abelian gauge theory in this regime. The strong-coupling expansion converges
  only for $\beta$ below a small threshold $\beta_0$ (monograph §6), and the physical trajectory
  $\beta(a)\to\infty$ leaves that region immediately.

Any proof must connect the two regimes through the crossover where neither perturbation theory nor the
strong-coupling expansion converges. The research program (`docs/03-research-program.md`) is organised
around that gap.

## 1.6 A quantitative restatement on the lattice

Let $m_L(\beta)$ denote the gap of the lattice Hamiltonian on a spatial lattice of $L^3$ sites in lattice
units, and $m(\beta)=\liminf_{L\to\infty}m_L(\beta)$. Two-loop asymptotic scaling for the Wilson action gives the
dimensionless function
$$a\Lambda_L(\beta)=\bigl(b_0g^2\bigr)^{-b_1/(2b_0^2)}\exp\Bigl(-\frac{1}{2b_0g^2}\Bigr),\qquad
b_0=\frac{11N}{48\pi^2},\quad b_1=\frac{34N^2}{3(16\pi^2)^2},\quad g^2=\frac{2N}{\beta}.$$
The **lattice mass gap conjecture** is the pair of statements

* **(G1)** $m(\beta)>0$ for every $\beta\in(0,\infty)$ (no massless phase at any coupling);
* **(G2)** $\liminf_{\beta\to\infty} m(\beta)/a\Lambda_L(\beta) > 0$ (the gap survives the continuum limit).

(G2) is the quantitative heart of the problem. (G1) alone would follow from (G2) plus the strong-coupling
result only if one also controlled the intermediate region; (G1) is also the lattice form of "confinement at
all couplings", and is open. The abelian theory $G=U(1)$ shows that neither statement can be proved by an
argument blind to the non-abelian structure: for $U(1)$ in four dimensions, (G1) is **false** at large
$\beta$ (Guth 1980; Fröhlich–Spencer 1982: a massless Coulomb phase exists).
