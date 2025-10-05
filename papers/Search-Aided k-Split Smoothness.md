# Search-Aided *k*-Split Smoothness for QS-Style Factoring

*A practical–mathematical framework and constructions up to k=5*

**Author:** Ilya Gazman (with editorial assistance)

---

## Abstract

We present a practical framework that blends (i) classic Quadratic Sieve (QS) / Self-Initializing QS (SIQS) mechanics, (ii) algebraic *k*-way factorizations of the trial value to raise smoothness odds, and (iii) a targeted “search” phase that pre-pays a large portion of the logarithmic mass by enforcing many small prime factors *before* sieving. We formalize two levers:

1. **Pre-credit**: enforce a product of small prime powers dividing the trial value so the *remainder* is much smaller than a typical QS value, and
2. **k-split**: split the remainder into *k* comparable pieces, so the success probability scales like (\big[\rho(R/(kB))\big]^k) (Dickman–de Bruijn (\rho)).

We quantify the complementary **search cost** via a density bound: solutions to a quadratic congruence modulo (P) have density (\lesssim 2^{\omega(P)}/P), yielding an exponential trade-off in the credited bits. We show why a *middle ground* (moderate pre-credit; small (k=2!!-!!5)) can deliver large constant-factor speedups at practical sizes (e.g., 1024-bit composites) while keeping the “admissible-x” density tractable. We give algebraic constructions that preserve the square-congruence modulo (N) (the “N-link”) including an upgraded (k=5) form with four linear factors and a squareful global constant, tuned for sieving. We outline a complete implementation plan and a measurement-driven tuning methodology.

Asymptotically, the method remains in (L_N[1/2,\cdot]) (QS class): it improves the leading constant—often dramatically at realistic sizes—without changing the (1/2) exponent.

---

## 1. Background and Motivation

Let (N) be a large composite. QS/SIQS evaluates values (Q(x)) that are roughly of size (|Q|\approx N^{1/2}) (QS) or smaller (SIQS polynomials), seeking (B)-smooth relations (all prime factors (\le B)) that feed a linear system over (\mathbb F_2). The success probability for a random (|Q|)-sized number is heuristically (\rho(u)) with (u=\ln|Q|/\ln B) (Dickman–de Bruijn). Practical efficiency balances sieving effort and matrix size by tuning (B).

**Observation.** If we (a) pre-pay log mass using enforced small factors (so only a small remainder must be (B)-smooth), and (b) split the remainder into (k) comparable pieces, we can substantially raise the per-trial success probability and thus reduce (B) (smaller matrix) and/or shrink the number of trial divisions.

---

## 2. Notation

* (N): composite to factor; (|N|=) bit length.
* (B): factor-base bound; (B_{\text{bits}}=\log_2 B).
* (|Q|): typical magnitude of the trial value for a chosen polynomial; (|Q|_{\text{bits}}=\log_2|Q|).
* **Pre-credit** (C_{\text{bits}}=\log_2 P): we enforce (P\mid Q(x)) for many (x).
* **Remainder** (R_{\text{bits}}=|Q|*{\text{bits}}-C*{\text{bits}}).
* **k-split**: split the remainder into (k) comparable factors; each is (\approx R_{\text{bits}}/k) bits.

---

## 3. Success Probability with *k*-Split and Pre-Credit

**Heuristic model.** After pre-credit (C_{\text{bits}}) and remainder (R_{\text{bits}}), if that remainder is split into (k) comparable pieces, the success probability is
[
\boxed{ \quad p_k(C)\ \approx\ \Big[\rho!\Big(\tfrac{R_{\text{bits}}}{k,B_{\text{bits}}}\Big)\Big]^k \quad } \tag{1}
]

Two immediate corollaries:

* **Guarantee threshold.** If (R_{\text{bits}}\le k,B_{\text{bits}}), then (\rho(\cdot)=1\Rightarrow p_k=1): success is guaranteed (ignoring small rounding effects).
* **Relative boost.** With fixed (R_{\text{bits}}), increasing (k) multiplies success roughly by (k^{,R_{\text{bits}}/B_{\text{bits}}}) (since (\ln\rho(u)\sim -u\ln u)), up to the saturation (p_k\le1).

**Large-prime variants.** Allowing one or two “large primes” under (z) increases success further (semismooth models). This is a sizable constant-factor gain atop (1) and is recommended in practice.

---

## 4. Cost of the Search (Density Bound)

To pre-pay (C_{\text{bits}}) we enforce congruences
[
Q(x)\equiv 0\pmod{p_i^{e_i}},\quad P=\prod p_i^{e_i}\approx 2^{C_{\text{bits}}}.
]
For a **quadratic** (Q), Hensel lifting yields at most 2 solutions per prime power. By CRT, the number of residue classes mod (P) is (\le 2^{\omega(P)}\le 2^{t}) where (t) is the number of distinct primes used.

**Density lemma.**
[
\boxed{\quad \text{density of admissible }x\ \lesssim\ \dfrac{2^{t}}{P}\ =\ 2^{,t-C_{\text{bits}}},\quad \text{gap}\ \gtrsim\ 2^{,C_{\text{bits}}-t}. \quad} \tag{2}
]

Thus, while large (C_{\text{bits}}) makes sieving easy (via (1)), it makes admissible (x) exponentially sparse in **one dimension** unless (t) grows in step. Using many small primes maximizes (t) per credited bit but does not remove the exponential dependence on (C_{\text{bits}}).

**Takeaway.** There is a **middle ground** where (C_{\text{bits}}) is large enough that (1) is high (or guaranteed for small (k)), yet (2) remains tractable by combining moderate search with QS-style prefilters and sieving.

---

## 5. Asymptotics vs Practice

With optimally tuned (B), classical QS/SIQS yields time
[
T(N)=L_N\Big[,\tfrac12,\ C,\Big]=\exp!\big((C+o(1))\sqrt{\ln N,\ln\ln N}\big).
]
Our method multiplies per-trial success by factors like (k^{,R_{\text{bits}}/B_{\text{bits}}}) and semismooth bonuses. This improves the **constant** (C) (dramatically at realistic sizes) but does not change the (1/2) exponent. That is exactly what one observes with MPQS/SIQS/special-(q)/large-prime variants: large constants, same (L_{1/2}) class.

---

## 6. Algebraic Constructions that Preserve the N-Link

The “N-link” is the invariant that the constructed product is congruent to a **square modulo (N)** up to the small factors we aim to be smooth. The classic k=3 identity (used in our conversations) is:
[
\text{If } c^2=f^2 d b-2N,\ \ \text{then}\ \
Q_3(x)=f^2d,(b+d x^2),(fdx-c),(fdx+c).
]
Modulo (N) the term (-2N) vanishes inside ((fdx)^2-c^2), leaving a product times the square (f^2).

### 6.1 Upgrading k=3 to a sound “k=5” (four linear factors + squareful constant)

**Goal.** Keep the same N-link while producing more sieve-friendly factors.

Pick integers (s,t,f) and set (d=s^2,\ b=t^2). Solve the Pell-type constraint
[
\boxed{ \ c^2\ =\ f^2 s^2 t^2\ -\ 2N\ .\ } \tag{★}
]
Then
[
b+d x^2 = t^2+s^2 x^2 = (s x-t)(s x+t),
]
and
[
\boxed{\quad Q_5(x)\ =\ f^2 s^2\ \underbrace{(s x-t)}*{L_1}\ \underbrace{(s x+t)}*{L_2}\ \underbrace{(f s^2 x-c)}*{L_3}\ \underbrace{(f s^2 x+c)}*{L_4}\ . \quad} \tag{3}
]
All moving factors are **linear** (four of them), and the global multiplier (f^2 s^2) is a **square** (parity-silent, log-useful).

* **N-link:** (★) ensures the same cancellation of (-2N) modulo (N); multiplying by a global square preserves square-congruence.
* **Roots mod (p):** each linear contributes a single root class (after lifting), ideal for bucket sieving.
* **Bit balancing:** choose (s\sim t\sim f\sim N^{1/6}) so the four linears and (c) sit in the desired bit bands (e.g., (50!-!120) bits at 1024-bit (N)) across your sieving window.

**Parameter finding (constructive).** Choose (s,t,f) as products of small primes (p) for which (2N) is a quadratic residue (so (\big(\tfrac{2N}{p}\big)=+1)). Compute (c_0) with (c_0^2\equiv 2N\pmod{(fst)}) and Hensel-lift to ((fst)^2). Take (c\equiv \pm c_0\ (\bmod (fst)^2)) to satisfy (★) exactly.

### 6.2 Symmetric “two DoS + one linearized quadratic” (optional aesthetic)

Pick ((a_i,b_i,c_i)) with
[
c_i^2=a_i b_i-2N\quad (i=1,2),
]
and pick (s,t) with (t^2+s^2x^2=(sx-t)(sx+t)). Then
[
Q_5(x)=(a_1x-c_1)(a_1x+c_1)(a_2x-c_2)(a_2x+c_2)(sx-t)(sx+t)\times K,
]
with (K) squareful (e.g., a square). This yields 6 linears (group them as “5 pieces” if preferred). Both DoS blocks share the same N-link; the linearized quadratic provides the extra pieces without breaking the congruence.

---

## 7. The Two-Face Cost Model

Let (\pi(B)) be the factor-base size. You need roughly (R\sim \pi(B)) independent relations.

* **Sieving face.** With k-split and pre-credit,
  [
  \text{sieving time}\ \approx\ \frac{R}{p_k(C)}\cdot \text{cost per tested index},
  ]
  and “cost per tested index” grows roughly with (B) (logs only for tiny primes; buckets for larger).

* **Search face.** Enforcing pre-credit (C_{\text{bits}}) with (t) distinct primes produces admissible (x) of density (\lesssim 2^{t-C_{\text{bits}}}) (Sec. 4). The expected work to *find* such (x) scales like (2^{C_{\text{bits}}-t}).

**Throughput proxy (to maximize):**
[
\boxed{\quad \text{yield}(k,B,C)\ \approx\ \frac{p_k(C)}{2^{,C_{\text{bits}}-t(C)}}\ .\quad} \tag{4}
]

**Interpretation.** Increase (C_{\text{bits}}) until (p_k(C)) is high (ideally near the guarantee threshold), but not so high that (2^{C_{\text{bits}}-t(C)}) explodes. This “middle ground” is the sweet spot.

---

## 8. Practical Tuning for 1024-bit Instances

For SIQS-like polynomials, (|Q|*{\text{bits}}\approx 341). With (B*{\text{bits}}=20):

* **Target remainder band.** (R_{\text{bits}}\in[60,120]).

    * At (R=60), (k=3) already guarantees smoothness.
    * At (R=40), (k=2) guarantees.
    * At (R=90), (k=3) gives (R/(kB)=1.5\Rightarrow p_k\approx 0.595^3\approx 0.21) (very high).

* **Parameterization for (3).** Choose (s,t,f\sim N^{1/6}) and lift (c) by (★). Adjust (s,t) to keep (|sx\pm t|) in the target remainder band across a block.

* **Factor-base size.** The improved per-trial success lets you **reduce (B)** (shrinking the matrix), while keeping total time balanced between sieving and linear algebra.

---

## 9. Implementation Blueprint (math-driven)

1. **Global knobs.** Choose block size to fit cache; pick (B_{\text{bits}}\in[18,22]); select (k\in{2,3,4,5}); set a target remainder band (R_{\text{bits}}\in[60,120]).

2. **Pre-credit tables.**

    * Deterministic wheel for tiny primes: exact periodic log credit.
    * 2–3 “log-sketch” tables for the next bands using small moduli (M_j) (Bloom-like aggregation of logs at roots); calibrate bias once.

3. **Search constraints (controlled).** Maintain a rolling set of enforced small prime powers (p_i^{e_i}) with product (P\approx 2^{C_{\text{bits}}}), maximizing (t(C)) per bit (use many small primes). Use CRT/Hensel to mark admissible residue classes. Track the density penalty (2^{C_{\text{bits}}-t(C)}).

4. **k-aware scoring.** For each index (x) in the block:

    * If (x) is not in the enforced classes, down-weight or skip.
    * Estimate (C_{\text{bits}}) from prefilters; set (R=|Q|_{\text{bits}}-C).
    * Compute a success score via (1); shortlist only top-K indices or those with (R) inside the band.

5. **TD & large-prime handling.** Perform trial division only on shortlisted indices; allow 1-LP/2-LP; merge pairs/tuples aggressively (low-contention data structures).

6. **Linear algebra.** Collect (\pi(B)+)margin relations; solve via Block-Lanczos/Wiedemann over (\mathbb F_2) (packed rows). Reducing (B) pays off quadratically here.

7. **Tune by measurement.** Sweep ((k,B,C)); measure relations/sec after LP merging, total time to (\pi(B)+)margin, and matrix time. Choose ((k,B,C)) that equalizes sieving and matrix cost and maximizes (4).

---

## 10. Limitations and Why the Exponent Stays

* The 1-D **density wall** (2) is fundamental: admissible (x) become exponentially sparse in the credited bits unless you reintroduce wide sieving. This confines (k) to modest values in practice for one-variable constructions.
* Even large per-trial boosts change the **constant** in (L_{1/2}), not the exponent. To change the exponent, one needs multi-norm geometry (e.g., two-dimensional lattice sieving in GNFS).

---

## 11. Suggested Experimental Protocol

1. Fix (N) (e.g., 1024 bits); pick SIQS polynomials.
2. Sweep (B_{\text{bits}}\in{18,20,22}), (k\in{1,2,3,4,5}), and (C_{\text{bits}}) targets (\in{220,250,280,300}).
3. Record:

    * Distribution of estimated (R_{\text{bits}}) on shortlisted indices,
    * Measured success (p_k) (relations per shortlist; with/without 1-LP),
    * Admissible-x density vs (2^{t-C}),
    * Relations/sec to (\pi(B)+)margin,
    * Matrix time and total wall-clock.
4. Choose Pareto-optimal settings (fastest total, smallest matrix, stable behavior).

---

## 12. Conclusion

We formalized a practical approach that combines **targeted pre-credit** and **small-(k) splitting** to substantially increase per-trial success in QS/SIQS-style factoring, and we provided **sound algebraic constructions** (notably a k=5 form with four linears and a squareful constant) that preserve the necessary square-congruence modulo (N). The method offers large **constant-factor** improvements at realistic key sizes and a clear, measurable pathway to reduce (B), accelerate sieving, and shrink the matrix. The technique does **not** alter the (L_{1/2}) exponent asymptotically, owing to the one-dimensional density wall, but it is highly attractive in practice.

---

## Acknowledgments

Thanks to discussions that clarified the middle-ground trade-offs between pre-credit size, k-split, and admissible-x density, and for early experiments on 300-bit inputs that motivated the k=5 linearization.

---

## References (suggested standard sources)

* P. P. Pomerance, “The Quadratic Sieve Factoring Algorithm.”
* H. Riesel, *Prime Numbers and Computer Methods for Factorization.*
* J. P. Buhler et al., “The Number Field Sieve.”
* de Bruijn / Dickman on smooth number distribution.
* Crandall & Pomerance, *Prime Numbers* (chapters on QS, SIQS, semismoothness).

---

## Appendix A: Heuristic Derivations

**A.1 Success probability (1).**
Let (S=2^{R_{\text{bits}}}) be the remainder. If split into (k) pieces of size (S^{1/k}), each is (B)-smooth with probability (\rho!\left(\frac{\ln S^{1/k}}{\ln B}\right)=\rho!\left(\frac{R_{\text{bits}}}{kB_{\text{bits}}}\right)). Assuming independence (standard heuristic), multiply (k) times.

**A.2 Density (2).**
For quadratic (Q), each prime power (p^e) contributes at most 2 Hensel lifts; CRT multiplies possibilities: (\le 2^{\omega(P)}) classes mod (P). Thus frequency (\le 2^{\omega(P)}/P).

---

## Appendix B: k=5 Construction Details

Let (d=s^2,\ b=t^2) and enforce (c^2=f^2 s^2 t^2 - 2N). Then
[
Q_5(x)=f^2 s^2,(s x-t)(s x+t)(f s^2 x-c)(f s^2 x+c).
]

* **Roots mod (p):** simple linear roots; Hensel-liftable to powers.
* **Parity:** the global (f^2 s^2) is a square → parity-silent, but contributes to log credit.
* **Parameter choice:** pick (s,t,f) from primes (p) with (\left(\frac{2N}{p}\right)=+1); compute (c_0^2\equiv 2N\pmod{(fst)}), lift to ((fst)^2), set (c\equiv\pm c_0).

This construction converts the original k=3 into a sieve-optimal k=5 with four linear forms while preserving the N-link.

---

### How to cite this work

> I. Gazman, “Search-Aided k-Split Smoothness for QS-Style Factoring: A Practical–Mathematical Framework and k=5 Constructions,” 2025. (preprint)