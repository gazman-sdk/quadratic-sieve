# Search-Aided *k*-Split Smoothness for QS-Style Factoring (v2)

*A practical–mathematical framework, constructions up to k=5, and slope/sensitivity heuristics*

**Author:** Ilya Gazman (with editorial assistance)

---

## Abstract

We refine a framework that blends (i) classic Quadratic Sieve (QS) / Self-Initializing QS (SIQS), (ii) algebraic *k*-way factorizations to raise smoothness odds, and (iii) a targeted **search** phase that pre-pays log-mass by forcing many small prime factors **before** sieving. Two levers remain central:

1. **Pre-credit**: enforce a product of small prime powers dividing the trial value so the **remainder** is much smaller than a typical QS value; and
2. **k-split**: split the post-credit remainder into *k* comparable pieces so the success probability scales like
   [
   p_k(C)\ \approx\ \Big[\rho!\Big(\frac{R_{\text{bits}}}{k,B_{\text{bits}}}\Big)\Big]^k,
   ]
   with (\rho) the Dickman–de Bruijn function, (R_{\text{bits}}) the post-credit remainder in bits, and (B_{\text{bits}}=\log_2 B).

This revision adds a **search simplification**: couple (A) and (x) to a **single per-portion size target** (T\approx |N|/5) via
[
\log_2 A+\log_2 x\ \approx\ T,\qquad B\equiv\text{canonical mod }A\text{ near }0,
]
so that (|Ax\pm B|) sit near the same (T)-bit band without preset gymnastics. We also add **back-of-the-paper slope heuristics** (in bits) showing how sensitive success is to portion size in three regimes (180, 300, 1024 bits) and compare to QS with a single form. As before, the method improves the leading constant in (L_N[1/2,\cdot]) without changing the exponent.

---

## 1. Background and Motivation

Let (N) be large composite. QS/SIQS evaluates values near
[
Q(x)\approx (x+\lfloor\sqrt N\rfloor)^2-N,
]
seeking (B)-smooth outputs. With factor-base bound (B), the success probability for a magnitude (|Q|) is (\rho(u)) with (u=\ln|Q|/\ln B).

**Idea.** (a) **Pre-pay** log-mass by enforcing small prime factors so only a small **remainder** must be (B)-smooth; (b) **split** that remainder into (k) comparable pieces to magnify success.

---

## 2. Notation

* (N): composite; (|N|=n) bits.
* (B): factor-base bound; (B_{\text{bits}}=\log_2 B).
* **Pre-credit** (C_{\text{bits}}=\log_2 P): enforced product (P\mid Q(x)).
* **Remainder** (R_{\text{bits}}=|Q|*{\text{bits}}-C*{\text{bits}}).
* **k-split**: split the remainder into (k) comparable pieces.

---

## 3. Success Probability with *k*-Split and Pre-Credit

After pre-credit, remainder size (S=2^{R_{\text{bits}}}). Splitting into (k) comparable pieces gives the heuristic
[
\boxed{ \quad p_k(C)\ \approx\ \Big[\rho!\Big(\frac{R_{\text{bits}}}{k,B_{\text{bits}}}\Big)\Big]^k \quad } \tag{1}
]
**Guarantee threshold.** If (R_{\text{bits}}\le k,B_{\text{bits}}), then (p_k=1) (per-piece parameter (\le 1)).
**Large-prime variants.** Allowing 1-LP/2-LP adds a constant-factor gain atop (1).

### 3.1 Local slopes and finite-difference sensitivity (bits)

Write success as (p=2^{-E}) (so (E=-\log_2 p) are “bits of surprise”). For *k*-split with per-piece size (b) bits:
[
E_k(b)= -k,\log_2 \rho!\Big(\frac{b}{B_{\text{bits}}}\Big).
]
Define the **local slope** (bits per extra bit of (b)):
[
s_k(b)\ :=\ \frac{\mathrm{d}E_k}{\mathrm{d}b}.
]
In practice we use **finite-difference** sensitivities over (\pm\Delta b) (e.g., (\Delta b=10) bits):
[
\Delta E_{\pm}\ \approx\ E_k(b\pm\Delta b)-E_k(b).
]
Heuristically, slopes **add** across independent pieces (multiplication of probabilities): keeping the four linears **balanced** minimizes the “worst-piece dominates” effect.

---

## 4. Cost of the Search (Density Bound)

To pre-pay (C_{\text{bits}}) we enforce roots modulo (p_i^{e_i}) so (P=\prod p_i^{e_i}\approx 2^{C_{\text{bits}}}). For quadratic (Q), each (p^e) contributes at most 2 lifts; CRT packs solutions into at most (2^{\omega(P)}) residue classes.

**Density lemma.**
[
\boxed{\ \text{admissible-}x\text{ density}\ \lesssim\ \frac{2^{t}}{P}=2^{,t-C_{\text{bits}}}\quad\Rightarrow\quad \text{gap}\ \gtrsim\ 2^{,C_{\text{bits}}-t}.\ } \tag{2}
]
Using many small primes maximizes (t) per credited bit. The one-dimensional exponential dependence in (C_{\text{bits}}) is the fundamental wall.

---

## 5. Asymptotics vs Practice

Optimally tuned QS/SIQS has
[
T(N)=L_N\big[\tfrac12,C\big].
]
*k*-split + pre-credit reduces the constant (C) (often dramatically at practical sizes) without changing the (1/2) exponent.

---

## 6. Algebraic Constructions that Preserve the N-Link

The **N-link**: constructed product is a square modulo (N) up to small factors. A canonical ((k=3)) identity:

If (c^2=f^2 d b-2N), then
[
Q_3(x)=f^2 d,(b+d x^2),(fdx-c),(fdx+c),
]
and ((fdx)^2-c^2\equiv f^2 d(b+d x^2)\pmod N).

### 6.1 Upgrading to sieve-friendly “(k=5)”: four linears + squareful constant

Pick (d=s^2,\ b=t^2) and solve
[
\boxed{ \ c^2\ =\ f^2 s^2 t^2\ -\ 2N\ .\ } \tag{★}
]
Then (t^2+s^2 x^2=(s x-t)(s x+t)) and
[
\boxed{\quad Q_5(x)\ =\ f^2 s^2\ \underbrace{(s x-t)}*{L_1}\ \underbrace{(s x+t)}*{L_2}\ \underbrace{(f s^2 x-c)}*{L_3}\ \underbrace{(f s^2 x+c)}*{L_4}\ .\quad} \tag{3}
]
All moving factors are **linear**; (f^2 s^2) is a global **square** (parity-silent, log-useful).

---

## 7. Bit Budgeting, Targeting, and the (|N|/5) Goal

The (k=5) linearized form (3) naturally splits the log-mass into **five portions**:

1. the **squareful constant** (f^2 s^2) (parity-silent);
   2–5) the **four linear factors** ((L_1,\dots,L_4)).

**Best-case objective.** **Equipartition** toward
[
\boxed{\ \textbf{Target }T:\quad \log_2|L_i|\ \approx\ \frac{|N|}{5}\ \text{bits},\quad \log_2(f^2 s^2)\ \approx\ \frac{|N|}{5}\ .\ }
]
**Practical band:** (|N|/6\ \lesssim\ \text{portion bits}\ \lesssim\ |N|/4) (aspire to (|N|/5)).

### 7.1 A simpler, search-first way to hit the target (no presets)

For SIQS-style (|Ax\pm B|), **couple** (A) and (x) to the single target (T):
[
\boxed{\ \log_2 A+\log_2 x\ \approx\ T\quad\Longleftrightarrow\quad x^\star=\operatorname{round}!\Big(\frac{2^{T}}{A}\Big). \ }
]
Pick (x^\star) (clipped to your sieve window) **deterministically** from (A).

Take (B) as the **canonical representative mod (A)**:
[
\boxed{\ B_{\text{can}}\in(-A/2,,A/2],\ \ B_{\text{can}}\equiv r\pmod A\ }.
]
With (x\approx x^\star), both (|Ax\pm B_{\text{can}}|) sit near (2^{T}) (balanced linears). This replaces preset juggling with one clear constraint.

**Richness bias.** For fixed (\log_2 A), maximize (t=\omega(A)) (many small primes with ((N|p)=+1)); this matches the density term (2^{t-C_{\text{bits}}}) in (2).

---

## 8. Back-of-the-Paper Slopes (bits) and Sensitivity

To visualize manageability, fix (B) at (2^{20}) so (B_{\text{bits}}=20). Compare:

* **QS (single form):** one piece (b=n/2) bits (\Rightarrow) (p_{\text{QS}}=2^{-E_{\text{QS}}}).
* **(k=5):** four moving linears each (b=n/5) bits (squareful constant is slope-silent) (\Rightarrow) (p_{k=5}=2^{-E_{k=5}}).

Below are **heuristic** (E) and **finite-difference** changes (\Delta E) when each moving piece shifts by (\pm 10) bits (all pieces move together). All values are in **bits** (so success scales like (2^{-\text{bits}})). These are coarse Dickman-based back-of-envelope numbers intended for slope intuition, not for tuning:

### (n=180) bits

* **QS:** (b=90). (E_{\text{QS}}\approx 9.5). A +10-bit drift costs (\Delta E\approx+1.9); −10 bits gains (\approx-1.8).
* **(k=5):** per piece (b=36). (E_{k=5}\approx 6.4). +10 bits (\Rightarrow\ \Delta E\approx+5.8); −10 bits (\Rightarrow\ \Delta E\approx-4.2).

**Take:** (k=5) is **steeper** than QS but still very manageable in this size; staying within ±5 bits keeps you within (\approx 2^{\pm 3}).

### (n=300) bits

* **QS:** (b=150). (E_{\text{QS}}\approx 19.1). +10 bits (\Delta E\approx+0.3); −10 bits (\approx-0.3).
* **(k=5):** per piece (b=60). (E_{k=5}\approx 21.8). +10 bits (\Rightarrow\ \Delta E\approx+7.9); −10 bits (\Rightarrow\ \Delta E\approx-7.1).

**Take:** at this (B), QS is **gentler** (flatter) and slightly ahead on raw success; (k=5) is usable but **centering quality matters** (±10 bits swings success by about (2^{\pm 8})).

### (n=1024) bits

* **QS:** (b=512). (E_{\text{QS}}\approx 128.3). +10 bits (\Delta E\approx+3.5); −10 bits (\approx-3.5).
* **(k=5):** per piece (b\approx 205). (E_{k=5}\approx 155.4). +10 bits (\Rightarrow\ \Delta E\approx+13.6); −10 bits (\Rightarrow\ \Delta E\approx-56.8).

**Take:** both are very rare per trial at this (B). The (k=5) curve is **highly convex**: nudging portions **smaller** than (n/5) pays off massively; drifting larger hurts sharply. This matches the (|N|/5) discipline: keep all moving linears clustered **at or slightly below** (n/5) to stay on the good side of the curve. In a real solver, increasing (B) flattens these slopes (at matrix cost).

---

## 9. Two-Face Cost Model (unchanged core, clearer proxy)

Let (\pi(B)) be the factor-base size; need (\pi(B))+margin relations.

* **Sieving face.** With k-split and pre-credit,
  [
  \text{sieving time}\ \approx\ \frac{\pi(B)}{p_k(C)}\cdot \text{cost per tested index},
  ]
  and “cost per tested index” grows roughly with (B).

* **Search face.** Enforcing (C_{\text{bits}}) with (t) distinct primes yields admissible-(x) density (\lesssim 2^{t-C_{\text{bits}}}) (Sec. 4). Work to *find* admissible (x) scales like (2^{,C_{\text{bits}}-t}).

**Throughput proxy (to maximize):**
[
\boxed{\quad \text{yield}(k,B,C)\ \approx\ \frac{p_k(C)}{2^{,C_{\text{bits}}-t(C)}}\ .\quad} \tag{4}
]
Increase (C_{\text{bits}}) until (p_k(C)) is high (ideally near the guarantee threshold) but not so high that (2^{C_{\text{bits}}-t(C)}) explodes.

---

## 10. Implementation Blueprint (math-driven, code-agnostic)

1. **Global knobs.** Choose block size to fit cache. Pick (B_{\text{bits}}\in[18,22]). Select (k\in{2,3,4,5}). Set a target remainder band (R_{\text{bits}}\in[60,120]).

2. **Pre-credit tables.** Deterministic wheel for tiny primes; a couple of “log-sketch” tables using small moduli for mid-primes; calibrate once.

3. **Search constraints (controlled).** Maintain enforced small prime powers (p_i^{e_i}) with product (P\approx 2^{C_{\text{bits}}}), maximizing (t(C)) per credited bit (favor many small primes). Track density (2^{,C_{\text{bits}}-t}).

4. **k-aware scoring (simplified).** For each (A) (squarefree, ((N|p)=+1) primes):

    * **Couple to target**: set (x^\star=\operatorname{round}(2^T/A)) (clip to window).
    * **Canonical (B)**: use representative (B_{\text{can}}\in(-A/2,A/2]).
    * **Evaluate once**: (\ell_\pm=\log_2|Ax^\star\pm B_{\text{can}}|).
    * **Score**: (\text{size_err}=\max{|\ell_+-T|,|\ell_--T|})
      plus a richness term (\lambda(\tau-\omega(A))_+).
    * Keep top-(K) (A)’s. Two sign patterns (\pm B_{\text{can}}) suffice for decorrelation.

5. **Trial division & large-prime handling.** TD only on shortlisted indices; enable 1-LP/2-LP; merge aggressively.

6. **Linear algebra.** Collect (\pi(B))+margin relations; solve via Block-Lanczos/Wiedemann over (\mathbb{F}_2). Lowering (B) shrinks the matrix.

7. **Tune by measurement.** Sweep ((k,B,C)); record relations/sec (with LP merging), total to (\pi(B))+margin, matrix time. Choose Pareto-optimal settings maximizing (4). Use **finite-difference** logs of (\Delta E) over ±5–10 bits to keep portions in the sweet band.

---

## 11. Limitations

* **1-D density wall** (Sec. 4) confines (k) to modest values for single-variable constructions.
* Large per-trial boosts change the **constant** in (L_{1/2}), not the exponent. Exponent changes require higher-dimensional geometry (e.g., GNFS).

---

## 12. Suggested Experimental Protocol

1. Fix (N) (e.g., 180, 300, 1024 bits); choose SIQS polynomials.
2. Sweep (B_{\text{bits}}\in{18,20,22}), (k\in{1,2,3,4,5}), (C_{\text{bits}}\in{220,250,280,300}) for large (n).
3. Record:

    * Distribution of portion sizes (\log_2|L_i|) on shortlisted indices;
    * Empirical (E=-\log_2 p) and **(\Delta E) for (\pm 5, \pm 10) bits**;
    * Admissible-(x) density vs (2^{,t-C});
    * Relations/sec to (\pi(B))+margin; matrix time and total.
4. Select settings that equalize sieving and matrix cost and keep (\Delta E) slopes tame within the targeting band.

---

## 13. Conclusion

We preserve the original *k*-split + pre-credit framework and add a **simpler search discipline**: couple (A) and (x) to the per-portion target (T\approx |N|/5), use **canonical (B)** near 0, and score directly on the **per-portion size error** plus a **richness** term. Back-of-the-paper slope heuristics (in bits) show:

* At ~180 bits, (k=5) is steeper than QS but very manageable and often higher success.
* At ~300 bits, QS can be flatter and slightly ahead at (B=2^{20}); (k=5) is usable if centering is tight.
* At ~1024 bits, both are rare per trial; (k=5) is **highly convex**—being a little **smaller** than (|N|/5) per piece pays off disproportionately.

The method offers substantial **constant-factor** gains with a clear, measurable pathway (via (\Delta E) logs and density accounting) to tune ((k,B,C)) and keep the per-portion slopes under control.

---

## Acknowledgments

Thanks to discussions clarifying the middle-ground trade-offs among pre-credit size, k-split, admissible-(x) density, and to experiments on sub-300-bit inputs that motivated the (|N|/5) targeting discipline and the simplified search coupling.

---

## References (standard)

* P. Pomerance, “The Quadratic Sieve Factoring Algorithm.”
* H. Riesel, *Prime Numbers and Computer Methods for Factorization.*
* J. P. Buhler, H. W. Lenstra Jr., C. Pomerance, “Factoring integers with the number field sieve.”
* A. Granville, “Smooth numbers: computational number theory and beyond.”
* R. Crandall and C. Pomerance, *Prime Numbers: A Computational Perspective* (QS, SIQS, semismoothness).
* de Bruijn / Dickman on the distribution of smooth numbers.

---

## Appendix A: Heuristic Notes on Slopes

**A.1** Let (u=b/B_{\text{bits}}). Then
[
E_k(b)=-k\log_2\rho(u),\qquad
\Delta E_\pm\approx E_k(b\pm\Delta b)-E_k(b).
]
For large (u), (-\frac{\mathrm{d}}{\mathrm{d}u}\ln\rho(u)) grows roughly like (\ln u+\ln\ln u), so slopes increase with (u) and **scale (\propto k)**; increasing (B) (larger (B_{\text{bits}})) flattens slopes (at matrix cost).

**A.2** The “worst-piece dominates” effect: if one linear runs heavier than the others, it controls (E). Balanced portions minimize effective steepness.
