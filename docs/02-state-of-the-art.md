# 2. State of the art

A map of what is rigorously known, what is known numerically, and what is physically believed.
Entries are labelled **[T]** theorem (rigorous, published), **[N]** numerical, **[P]** physical
argument, **[C]** claimed but not generally accepted. References are listed at the end; those we
have not yet re-read in full during this project are marked (verify) in the bibliography and must be
checked before being relied on in a proof.

## 2.1 Rigorous results on lattice gauge theory

| Statement | Status | Where |
|---|---|---|
| Wilson's lattice measure is a well-defined probability measure for any compact $G$, any $\beta$, any finite lattice; gauge invariance; Elitzur's theorem (local gauge symmetry cannot break) | [T] | Wilson 1974; Elitzur 1975; monograph §2 |
| Reflection positivity of the Wilson action for site reflections (and link reflections) | [T] | Osterwalder–Seiler 1978; monograph §3 |
| Self-adjoint, strictly positive transfer matrix; Hamiltonian; unique vacuum on a finite lattice | [T] | Lüscher 1977; Osterwalder–Seiler 1978; monograph §4 |
| Strong-coupling (small $\beta$) convergent cluster expansion: exponential clustering uniform in volume, area law for Wilson loops, mass gap, analyticity in $\beta$, existence of the infinite-volume limit | [T] | Osterwalder–Seiler 1978; Seiler 1982; monograph §6; in this project: Part II (infinite volume, uniqueness, Markov property, OS reconstruction, gap; $\beta<\operatorname{artanh}\frac1{6(d-1)}$), Part III (area law), Part IV (gap between $4\log(1/u(\beta))-\log(1+\delta(\beta))-O(\beta^5)$ and $4\log(1/u(\beta))+D(\beta)$, explicit, for all $\beta<\beta_0(d)$, $d\ge3$; to first order $m=4\log(1/u)-\kappa(\beta)+O(\beta^2)$, with $\kappa=O(\beta^2)$ for $SU(2)$ and $\kappa=\beta/6+O(\beta^2)$ for $SU(3)$; to second order
$m=4\log(1/u)-\log\Upsilon(\beta)+O(\beta^3)$, $=4\log(1/u)+\beta^2/8+O(\beta^3)$ for $SU(2)$; to third order nothing new,
$m=4\log(1/u)-\log\Upsilon(\beta)+O(\beta^4)$; at fourth order the upper bound $m\le4\log(1/u)-\log\Upsilon-(26d-70)u^4+o(\beta^4)$ for $SU(N)$,
$N\le5$, and the conjecture of equality; the lower bound at that order, $m\ge4\log(1/u)-\log\Upsilon-(\delta'_4+o(1))\beta^4$, with $\delta'_4$ explicit, about
$166$ for $SU(2)$ in $d=4$), and with the moves by characters the fourth-order term $m=4\log(1/u)-\log\Upsilon-(26d-70)u^4+o(\beta^4)$ for $SU(N)$, every $N$ (for $N\ge6$ via $\Upsilon_+\ge\Upsilon_-$ up to $O(\beta^5)$, Part IV Lemma 15.6),
and for every self-conjugate $\rho$ (and $-8u^4$ for every $\rho$ in $d=3$), refereed; the two sectors of charge conjugation separately (Part IV §15): each to third order, the even one to fourth, the $SU(3)$ splitting $\beta/3$; the odd sector at fourth order ($U(1)$, $d\ge4$) is a conjecture, with the hop/turn cancellation proved as an operator identity; closed polymers have 6 or $\ge8$ plaquettes, so every lower bound loses only $O(\beta^6)$, and with the junctions of three front-pairs of lowest order counted exactly the $SU(3)$ third-order term is resolved for $\beta\le1.7\cdot10^{-6}$ ($d=4$) and $\beta\le1.8\cdot10^{-5}$ ($d=3$) (Part IV §16, numerical evaluation of explicit constants); glueball mass $-4\log\beta+r(\beta)$, $r$ analytic: Schor 1983–84 (reported in Faria da Veiga–O'Carroll 2025; originals: verify); strong-coupling series of the glueball masses to
$u^8$ for $SU(2)$, $SU(3)$: Münster 1981, Seo 1982, as quoted by Langelage–Münster–Philipsen, PoS LAT2007:201 and JHEP 0807:036 (read; originals not) |
| Area law for $U(N)$ lattice Yang–Mills for $\beta\lesssim_d N$ (normalisation of the monograph; uniformly in $N$ in the normalisation $\exp(2\beta\,\mathrm{Re}\Tr U_p)$), via the master loop equation | [T] | Cao–Nissim–Sheffield, arXiv:2505.16585 (2025), Thm 1.2 (statement checked by a referee of Part III) |
| Area law in the 't Hooft regime $\beta\lesssim_d N^2$ (monograph normalisation) for $U(N)$, $SU(N)$, $SO(2N)$, by verifying the Durhuus–Fröhlich $\sigma$-model mass-gap criterion with Bakry–Émery methods | [T] | Cao–Nissim–Sheffield, arXiv:2509.04688 (2025) (threshold $\beta<1/(8(d-1))$ for the action $N\beta\,\mathrm{Re}\Tr$; checked by a referee of Part III) |
| Log-Sobolev inequality / exponential decay at strong coupling via Bakry–Émery, using the positive Ricci curvature of $SU(N)$, $SO(N)$ (a non-abelian input, though still a strong-coupling result) | [T] | Shen–Zhu–Zhu, Commun. Math. Phys. 400 (2023); arXiv:2204.12737 (verify the exact threshold) |
| $U(1)$ in $d=4$: a massless (Coulomb) phase exists at large $\beta$; Wilson loops obey a perimeter law | [T] | Guth 1980; Fröhlich–Spencer 1982 |
| $U(1)$ in $d=3$: confinement and mass gap for **all** $\beta$ (Polyakov's monopole mechanism made rigorous) | [T] | Göpfert–Mack 1982 |
| $\mathbb{Z}_2$ gauge theory in $d=4$ is dual to the Ising model; confinement–deconfinement transition | [T] | Wegner 1971; Balian–Drouffe–Itzykson 1975 |
| Finite-temperature deconfinement at large $\beta$ with few time slices ($SU(N)$) | [T] | Borgs–Seiler 1983 |
| Inequalities relating Wilson loops to vortex/electric-flux free energies ($SU(2)$) | [T] | Tomboulis–Yaffe 1985 |
| Two-dimensional Yang–Mills: exactly solvable; continuum limit; Wilson loop expectations for all $\beta$ | [T] | Migdal 1975; Driver 1989; Gross–King–Sengupta 1989; Sengupta 1997; Lévy 2010 |
| Large-$N$ strong-coupling: $SO(N)$, $SU(N)$ Wilson loops at strong coupling solved via string trajectories; $1/N$ expansion | [T] | Chatterjee 2019; Chatterjee–Jafarov 2016; Jafarov 2016; Basu–Ganguly 2018 |
| A criterion for the area law in terms of center-symmetry breaking (any $d$, $G$ with non-trivial center) | [T] | Chatterjee 2021 |
| Leading term of the free energy of 4-d $SU(N)$ lattice theory as $\beta\to\infty$ | [T] | Chatterjee 2016 |
| Ultraviolet stability of 4-d lattice Yang–Mills in finite volume: bounds on the partition function uniform in the lattice spacing along the asymptotic-freedom trajectory | [T] | Balaban 1984–1989 (series of ~15 papers); Magnen–Rivasseau–Sénéor 1993 |
| Confinement for all $\beta$ in 4-d $SU(2)$ | [C] | Tomboulis 2007 (gaps in the interpolation argument; not accepted) |
| Existence of the infinite-volume limit of 4-d $SU(N)$ lattice theory for all $\beta$ (convergence as $L\to\infty$, not just subsequential) | open | listed as open in Chatterjee's survey (verify) |
| Mass gap uniform in volume for 4-d $SU(N)$ at **large** $\beta$ | open | — |

## 2.2 Rigorous results in the continuum

| Statement | Status | Where |
|---|---|---|
| Osterwalder–Schrader reconstruction (Euclidean axioms ⇒ Wightman theory) | [T] | Osterwalder–Schrader 1973, 1975 |
| Constructive field theory in $d=2,3$: $P(\phi)_2$, $\phi^4_3$, Yukawa$_2$, Gross–Neveu$_2$; existence, OS axioms, mass gap at weak coupling | [T] | Glimm–Jaffe and many others, 1968–1985 |
| Abelian Higgs model in $d=2,3$: construction with OS axioms | [T] | Brydges–Fröhlich–Seiler 1979–1981 |
| 2-d Yang–Mills measure as a stochastic process on loops; Langevin dynamics for the 2-d Yang–Mills measure (regularity structures) | [T] | Lévy 2010; Chandra–Chevyrev–Hairer–Shen 2022 |
| 3-d Yang–Mills(–Higgs): local-in-time solutions of the stochastic quantisation equation; a state space for the 3-d measure | [T] | Chandra–Chevyrev–Hairer–Shen 2022; Cao–Chatterjee 2023 |
| Existence of the 3-d Yang–Mills measure in the continuum (finite volume) | open (parts of the UV control exist: Balaban's 3-d papers) | — |
| 4-d Yang–Mills: existence of the continuum measure, even in finite volume | open | — |
| Mass gap in 4-d | open | — |

## 2.3 Known no-go facts and necessary features of any proof

1. **Abelian obstruction.** $U(1)_4$ has a massless phase at weak coupling [T]. Hence no argument that uses
   only compactness of $G$, reflection positivity, and the structure of the Wilson action can prove (G1):
   the proof must use non-abelian features (non-trivial commutators, asymptotic freedom, the centre $\mathbb{Z}_N$,
   non-abelian "charge screening").
2. **Perturbation theory is blind to the gap.** All orders of lattice or continuum perturbation theory are
   massless. The gap is of the form $\Lambda\propto\exp(-1/(2b_0g^2))$, beyond all orders.
3. **Large volume is forced.** On a torus of size $L\ll 1/\Lambda$ the theory is weakly coupled and the gap
   is $O(1/L)$ (Lüscher 1983; Lüscher–Münster 1984, the "femto-universe"); the physical gap only appears for
   $L\gtrsim 1/\Lambda$, i.e. one cannot avoid the strong-coupling regime by working in small volumes.
4. **Two-dimensional analogue.** The $O(N)$ non-linear sigma model in $d=2$ is asymptotically free for $N\ge3$
   and is believed to have a mass gap for all couplings; this is also unproven, while for $N=2$ (XY model)
   the massless Kosterlitz–Thouless phase is a theorem (Fröhlich–Spencer 1981). The pair
   ($O(2)_2$ massless proven, $O(N\ge3)_2$ gap open) is the exact analogue of ($U(1)_4$ massless proven,
   $SU(N)_4$ gap open). A technique that proves the $O(3)$ gap would be a serious candidate for Yang–Mills.
   (Patrascioiu and Seiler have argued the opposite conjecture for $O(3)$; the rigorous status is open.)
5. **Bulk transitions are not obstacles.** For $SU(N)$ with $N\ge5$ (and for mixed fundamental–adjoint actions)
   the Wilson action has first-order *bulk* transitions as a function of $\beta$ [N]; both phases confine. The
   correct conjecture is (G1) for all $\beta$, not "no phase transition".

## 2.4 Numerical knowledge (what the answer looks like)

All of the following are [N], from lattice Monte Carlo over four decades.

* **Crossover, not transition** (for $SU(2)$, $SU(3)$): the plaquette interpolates smoothly between the
  strong-coupling series $\langle P\rangle=u(\beta)+16u^5u'(\beta)+\dots$ ($SU(2)$, $d=4$) and the weak-coupling series
  $1-(N^2-1)/(4\beta)+\dots$, with the crossover at $\beta\approx2.2$ ($SU(2)$) and $\beta\approx5.6$ ($SU(3)$).
  Our scans (`results/`) reproduce this.
* **String tension and glueballs.** The lightest state is the scalar glueball $0^{++}$ with
  $m_{0^{++}}/\sqrt\sigma\approx3.4$–$3.8$ for $SU(2)$ and $SU(3)$ (Teper 1998; Morningstar–Peardon 1999;
  Lucini–Teper–Wenger 2004; Athenodorou–Teper 2020), followed by $2^{++}$ at $\approx 1.4\,m_{0^{++}}$. In
  physical units $m_{0^{++}}\approx1.7$ GeV for $SU(3)$.
* **Asymptotic scaling.** $a\sqrt\sigma$ and $a\,m_{0^{++}}$ track (two-loop) asymptotic scaling once
  $\beta\gtrsim2.4$ ($SU(2)$) / $\beta\gtrsim6.0$ ($SU(3)$) up to known lattice artefacts; ratios of masses
  are constant in this region.
* **Finite-temperature deconfinement.** For $N_t$ time slices there is a transition (second order for $SU(2)$,
  first order for $SU(3)$) at $\beta_c(N_t)$ with $T_c/\sqrt\sigma\approx0.69$ ($SU(2)$), $0.63$ ($SU(3)$).
  Zero-temperature statements require $N_t\gg N_s$ or $N_t\gtrsim$ several correlation lengths.
* **Large $N$.** Glueball masses and $\sqrt\sigma$ have smooth $1/N^2$ expansions; the $N=\infty$ theory is
  confining with a gap.

These numbers tell us what the proof must reproduce: a gap of order $3.5\sqrt\sigma$ that persists as
$\beta\to\infty$ along the asymptotic-scaling curve.

## 2.5 Where the known techniques stop

| Technique | Regime it controls | Why it stops |
|---|---|---|
| Strong-coupling cluster expansion | $\beta<\beta_0\ll1$ | the expansion parameter $u(\beta)\to1$ as $\beta\to\infty$; no convergence beyond a finite radius |
| Dobrushin uniqueness | $\beta<\beta_D$ (even smaller) | same |
| Weak-coupling (Balaban) RG | $a\to0$ at fixed physical volume, small running coupling | at scale $1/\Lambda$ the running coupling is $O(1)$; small-field/large-field decomposition no longer effective |
| Stochastic quantisation / regularity structures | $d=2$ fully, $d=3$ locally in time | $d=4$ is critical for the Yang–Mills SPDE; theory not available |
| Exact solution / character expansion | $d=2$ | needs a topological (area-preserving) theory |
| Duality + Debye screening (Göpfert–Mack) | abelian, $d=3$ | no exact duality for non-abelian groups |
| Large $N$ | strong coupling only (rigorously) | planar limit at weak coupling not controlled |

## 2.6 Bibliography (verify = re-read before citing in a proof)

* Jaffe, A., Witten, E.: *Quantum Yang–Mills theory*. Clay Mathematics Institute Millennium Prize problem description (2000).
* Douglas, M.: *Report on the status of the Yang–Mills Millennium Prize problem* (2004). (verify)
* Wilson, K.: Confinement of quarks. Phys. Rev. D **10**, 2445 (1974).
* Osterwalder, K., Seiler, E.: Gauge field theories on a lattice. Ann. Phys. **110**, 440–471 (1978).
* Seiler, E.: *Gauge Theories as a Problem of Constructive Quantum Field Theory and Statistical Mechanics*. Lecture Notes in Physics **159**, Springer (1982).
* Lüscher, M.: Construction of a selfadjoint, strictly positive transfer matrix for Euclidean lattice gauge theories. Commun. Math. Phys. **54**, 283–292 (1977).
* Osterwalder, K., Schrader, R.: Axioms for Euclidean Green's functions I, II. Commun. Math. Phys. **31**, 83 (1973); **42**, 281 (1975).
* Streater, R., Wightman, A.: *PCT, Spin and Statistics, and All That* (1964).
* Glimm, J., Jaffe, A.: *Quantum Physics: A Functional Integral Point of View*, 2nd ed., Springer (1987).
* Elitzur, S.: Impossibility of spontaneously breaking local symmetries. Phys. Rev. D **12**, 3978 (1975).
* Kogut, J., Susskind, L.: Hamiltonian formulation of Wilson's lattice gauge theories. Phys. Rev. D **11**, 395 (1975).
* Guth, A.: Existence proof of a nonconfining phase in four-dimensional U(1) lattice gauge theory. Phys. Rev. D **21**, 2291 (1980).
* Fröhlich, J., Spencer, T.: Massless phases and symmetry restoration in abelian gauge theories and spin systems. Commun. Math. Phys. **83**, 411–454 (1982).
* Fröhlich, J., Spencer, T.: The Kosterlitz–Thouless transition in two-dimensional abelian spin systems and the Coulomb gas. Commun. Math. Phys. **81**, 527–602 (1981).
* Göpfert, M., Mack, G.: Proof of confinement of static quarks in 3-dimensional U(1) lattice gauge theory for all values of the coupling constant. Commun. Math. Phys. **82**, 545–606 (1982).
* Borgs, C., Seiler, E.: Lattice Yang–Mills theory at nonzero temperature and the confinement problem. Commun. Math. Phys. **91**, 329–380 (1983).
* Tomboulis, E., Yaffe, L.: Finite temperature SU(2) lattice gauge theory. Commun. Math. Phys. **100**, 313–341 (1985). (verify)
* Mack, G., Petkova, V.: Comparison of lattice gauge theories with gauge groups $\mathbb{Z}_2$ and SU(2). Ann. Phys. **123**, 442 (1979). (verify)
* Tomboulis, E.: Confinement for all values of the coupling in four-dimensional SU(2) gauge theory. arXiv:0707.2179 (2007). [C]
* Balaban, T.: Renormalization group approach to lattice gauge field theories I, II. Commun. Math. Phys. **109**, 249–301 (1987); **116**, 1–22 (1988); Convergent renormalization expansions for lattice gauge theories. Commun. Math. Phys. **119**, 243–285 (1988); Large field renormalization I, II. Commun. Math. Phys. **122**, 175–202, 355–392 (1989); and earlier papers 1984–1985. (verify exact list)
* Magnen, J., Rivasseau, V., Sénéor, R.: Construction of YM$_4$ with an infrared cutoff. Commun. Math. Phys. **155**, 325–383 (1993).
* Lüscher, M.: Some analytic results concerning the mass spectrum of Yang–Mills gauge theories on a torus. Nucl. Phys. B **219**, 233 (1983). Lüscher, M., Münster, G.: Nucl. Phys. B **232**, 445 (1984).
* Brydges, D., Fröhlich, J., Seiler, E.: On the construction of quantized gauge fields I–III. Ann. Phys. **121**, 227 (1979); Commun. Math. Phys. **71**, 159 (1980); **79**, 353 (1981).
* Driver, B.: YM$_2$: continuum expectations, lattice convergence, and lassos. Commun. Math. Phys. **123**, 575 (1989). Gross, L., King, C., Sengupta, A.: Two-dimensional Yang–Mills theory via stochastic differential equations. Ann. Phys. **194**, 65 (1989). Sengupta, A.: Gauge theory on compact surfaces. Mem. AMS **126** (1997). Lévy, T.: Two-dimensional Markovian holonomy fields. Astérisque **329** (2010).
* Chatterjee, S.: Yang–Mills for probabilists. In: *Probability and Analysis in Interacting Physical Systems*, Springer (2019). arXiv:1803.01950.
* Chatterjee, S.: Rigorous solution of strongly coupled SO(N) lattice gauge theory in the large N limit. Commun. Math. Phys. **366**, 203–268 (2019). Chatterjee, S., Jafarov, J.: The $1/N$ expansion for SO(N) lattice gauge theory at strong coupling. arXiv:1604.04777. Jafarov, J.: Wilson loop expectations in SU(N) lattice gauge theory. arXiv:1610.03821. Basu, R., Ganguly, S.: SO(N) lattice gauge theory, planar and beyond. CPAM (2018).
* Chatterjee, S.: A probabilistic mechanism for quark confinement. Commun. Math. Phys. **385**, 1007–1039 (2021).
* Chatterjee, S.: The leading term of the Yang–Mills free energy. J. Funct. Anal. **271**, 2944–3005 (2016).
* Cao, S., Chatterjee, S.: A state space for 3D Euclidean Yang–Mills theories. Commun. Math. Phys. (2023). arXiv:2111.12813.
* Chandra, A., Chevyrev, I., Hairer, M., Shen, H.: Langevin dynamic for the 2D Yang–Mills measure. Publ. Math. IHÉS **136**, 1–147 (2022); Stochastic quantisation of Yang–Mills–Higgs in 3D. arXiv:2201.03487 (2022).
* Cao, S., Nissim, R., Sheffield, S.: Expanded regimes of area law for lattice Yang–Mills theories. arXiv:2505.16585 (2025); Dynamical approach to area law for lattice Yang–Mills. arXiv:2509.04688 (2025).
* Shen, H., Zhu, R., Zhu, X.: A stochastic analysis approach to lattice Yang–Mills at strong coupling. Commun. Math. Phys. **400**, 805–851 (2023). arXiv:2204.12737. (verify)
* Fernández, R., Procacci, A.: Cluster expansion for abstract polymer models. New bounds from an old approach. Commun. Math. Phys. **274**, 123–140 (2007). arXiv:math-ph/0605041.
* Kotecký, R., Preiss, D.: Cluster expansion for abstract polymer models. Commun. Math. Phys. **103**, 491–498 (1986).
* Dobrushin, R.: The description of a random field by means of conditional probabilities and conditions of its regularity. Theory Probab. Appl. **13**, 197–224 (1968); Prescribing a system of random variables by conditional distributions. Theory Probab. Appl. **15**, 458–486 (1970).
* Friedli, S., Velenik, Y.: *Statistical Mechanics of Lattice Systems: A Concrete Mathematical Introduction*. Cambridge (2017). (cluster expansion, Dobrushin uniqueness, Peierls, etc.)
* Simon, B.: *The Statistical Mechanics of Lattice Gases*, Vol. I. Princeton (1993).
* Creutz, M.: Monte Carlo study of quantized SU(2) gauge theory. Phys. Rev. D **21**, 2308 (1980); *Quarks, Gluons and Lattices*, Cambridge (1983).
* Cabibbo, N., Marinari, E.: A new method for updating SU(N) matrices in computer simulations of gauge theories. Phys. Lett. B **119**, 387 (1982). Kennedy, A., Pendleton, B.: Improved heatbath method for Monte Carlo SU(2) gauge theory. Phys. Lett. B **156**, 393 (1985).
* Teper, M.: Glueball masses and other physical properties of SU(N) gauge theories in D=3+1: a review of lattice results for theorists. arXiv:hep-th/9812187 (1998). Morningstar, C., Peardon, M.: The glueball spectrum from an anisotropic lattice study. Phys. Rev. D **60**, 034509 (1999). Lucini, B., Teper, M., Wenger, U.: Glueballs and k-strings in SU(N) gauge theories. JHEP **06**, 012 (2004). Athenodorou, A., Teper, M.: The glueball spectrum of SU(3) gauge theory in 3+1 dimensions. JHEP **11**, 172 (2020).
* Kupiainen, A.: On the $1/n$ expansion. Commun. Math. Phys. **73**, 273–294 (1980). (verify)
* Patrascioiu, A., Seiler, E.: various papers 1990s–2000s questioning the $O(3)$ mass gap. (verify)
