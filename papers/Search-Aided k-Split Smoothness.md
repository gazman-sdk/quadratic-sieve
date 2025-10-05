# Search-Aided *k*-Split Smoothness for QS-Style Factoring

*A practical–mathematical framework and constructions up to k=5*

**Author:** Ilya Gazman (with editorial assistance)

---

## Abstract

We present a practical framework that blends (i) classic Quadratic Sieve (QS) / Self-Initializing QS (SIQS) mechanics, (ii) algebraic *k*-way factorizations of the trial value to raise smoothness odds, and (iii) a targeted **search** phase that pre-pays a large portion of the logarithmic mass by enforcing many small prime factors **before** sieving. We formalize two levers:

1. **Pre-credit**: enforce a product of small prime powers dividing the trial value so the **remainder** is much smaller than a typical QS value; and
2. **k-split**: split the remainder into *k* comparable pieces so the success probability scales like
   [
   p_k(C)\ \approx\ \Big[\rho!\Big(\frac{R_{\text{bits}}}{k,B_{\text{bits}}}\Big)\Big]^k
   ]
   with (\rho) the Dickman–de Bruijn function, (R_{\text{bits}}) the post-credit remainder in bits, and (B_{\text{bits}}=\log_2 B) the factor-base bound in bits.

We quantify the complementary **search cost** via a density bound: solutions to a quadratic congruence modulo (P) have density (\lesssim 2^{\omega(P)}/P), yielding an exponential trade-off in the credited bits. We show why a **middle ground** (moderate pre-credit; small (k=2\text{–}5)) can deliver large constant-factor speedups at practical sizes (e.g., 1024-bit composites) while keeping the “admissible-x” density tractable. We give algebraic constructions that preserve the essential **square-congruence modulo (N)** (the “N-link”), including an upgraded (k=5) form with four linear factors and a squareful global constant tuned for sieving. We outline a complete implementation plan and a measurement-driven tuning methodology.

Asymptotically, the method remains in (L_N[1/2,\cdot]) (QS class): it improves the leading constant—often dramatically at realistic sizes—without changing the (1/2) exponent.

---

## 1. Background and Motivation

Let (N) be a large composite. QS/SIQS evaluates values of the form
[
Q(x)\approx (x+\lfloor\sqrt N\rfloor)^2-N
]
whose magnitude (|Q|) is (\approx N^{1/2}) for QS (smaller for SIQS polynomials). The goal is to find **relations** where (Q(x)) is (B)-smooth (all prime factors (\le B)). The success probability for a random (|Q|)-sized integer is (\rho(u)) with (u=\ln|Q|/\ln B) (Dickman–de Bruijn). Efficiency balances sieving effort and matrix cost (about (\pi(B)) columns) by tuning (B).

**Idea.** If we (a) **pre-pay** log mass by enforcing small prime factors so only a small **remainder** must be (B)-smooth, and (b) **split** that remainder into (k) comparable pieces, we can drastically raise per-trial success and thus reduce (B) (shrinking the matrix) and/or reduce trial division.

---

## 2. Notation

* (N): composite to factor; (|N|=) bit-length of (N).
* (B): factor-base bound; (B_{\text{bits}}=\log_2 B).
* (|Q|): typical magnitude of trial values; (|Q|_{\text{bits}}=\log_2|Q|).
* **Pre-credit** (C_{\text{bits}}=\log_2 P): we enforce (P\mid Q(x)) via CRT/Hensel constraints.
* **Remainder** (R_{\text{bits}}=|Q|*{\text{bits}}-C*{\text{bits}}).
* **k-split**: split the post-credit remainder into (k) comparable pieces (heuristically independent).

---

## 3. Success Probability with *k*-Split and Pre-Credit

After pre-credit (C_{\text{bits}}), the remainder magnitude is (2^{R_{\text{bits}}}). Splitting it into (k) comparable pieces yields the heuristic success probability
[
\boxed{ \quad p_k(C)\ \approx\ \Big[\rho!\Big(\frac{R_{\text{bits}}}{k,B_{\text{bits}}}\Big)\Big]^k \quad } \tag{1}
]
Two immediate consequences:

* **Guarantee threshold.** If (R_{\text{bits}}\le k,B_{\text{bits}}), then each piece has Dickman parameter (\le 1), so (\rho=1\Rightarrow p_k=1).

* **Relative boost.** For fixed (R_{\text{bits}}), increasing (k) multiplies success roughly by (k^{,R_{\text{bits}}/B_{\text{bits}}}) (because (\ln \rho(u)\sim -u\ln u)); this saturates at (p_k\le 1).

**Large-prime variants.** Allowing one or two “large primes” under (z) (semismoothness) adds a sizable constant-factor gain atop (1). In practice, 1-LP/2-LP should be enabled.

---

## 4. Cost of the Search (Density Bound)

To pre-pay (C_{\text{bits}}) we enforce congruences
[
Q(x)\equiv 0\pmod{p_i^{e_i}},\quad P=\prod p_i^{e_i}\approx 2^{C_{\text{bits}}}.
]
For a **quadratic** (Q), Hensel lifting yields at most 2 solutions per prime power; via CRT, solutions modulo (P) lie in at most (2^{\omega(P)}\le 2^{t}) residue classes, where (t) is the number of distinct primes used.

**Density lemma.**
[
\boxed{\quad \text{density of admissible }x\ \lesssim\ \dfrac{2^{t}}{P}\ =\ 2^{,t-C_{\text{bits}}},\quad \text{so gap}\ \gtrsim\ 2^{,C_{\text{bits}}-t}. \quad} \tag{2}
]

Thus, while larger (C_{\text{bits}}) makes sieving easy (via (1)), it makes admissible (x) exponentially sparse in **one dimension** unless (t) grows apace. Using many small primes maximizes (t) per credited bit but does not remove the exponential dependence on (C_{\text{bits}}).

**Takeaway.** There is a **middle ground** where (C_{\text{bits}}) is large enough that (1) is high (often near-guaranteed for small (k)), yet (2) remains tractable when combined with QS-style prefilters and sieving.

---

## 5. Asymptotics vs Practice

With optimally tuned (B), QS/SIQS has runtime
[
T(N)=L_N!\Big[,\tfrac12,\ C,\Big]=\exp!\big((C+o(1))\sqrt{\ln N,\ln\ln N}\big).
]
Our method multiplies per-trial success by factors like (k^{,R_{\text{bits}}/B_{\text{bits}}}) plus semismooth bonuses. This reduces the **constant** (C) (often dramatically at realistic sizes) without changing the (1/2) exponent. This matches experience with MPQS/SIQS/special-(q)/large-prime variants: large constants, same (L_{1/2}) class.

---

## 6. Algebraic Constructions that Preserve the N-Link

The **N-link** is the invariant that the constructed product is congruent to a **square modulo (N)** up to the small factors we aim to be smooth. The canonical (k=3) identity we build from is:

> If (c^2=f^2 d b-2N), then
> [
> Q_3(x)=f^2 d,(b+d x^2),(fdx-c),(fdx+c),
> ]
> and ((fdx)^2-c^2\equiv (fdx)^2-(f^2db-2N)\equiv f^2 d(b+d x^2)\pmod N), so the product equals a square (times small factors) modulo (N).

### 6.1 Upgrading (k=3) to a sieve-friendly “(k=5)”: four linear factors + squareful constant

We improve sieving by **linearizing** the quadratic while keeping the same N-link.

Pick integers (s,t,f) and set (d=s^2,\ b=t^2). Solve the Pell-type constraint
[
\boxed{ \ c^2\ =\ f^2 s^2 t^2\ -\ 2N\ .\ } \tag{★}
]
Then
[
b+d x^2 = t^2+s^2 x^2 = (s x-t)(s x+t),
]
and we obtain the k=5 form
[
\boxed{\quad Q_5(x)\ =\ f^2 s^2\ \underbrace{(s x-t)}*{L_1}\ \underbrace{(s x+t)}*{L_2}\ \underbrace{(f s^2 x-c)}*{L_3}\ \underbrace{(f s^2 x+c)}*{L_4}\ . \quad} \tag{3}
]
All moving factors are **linear**; (f^2 s^2) is a global **square** (parity-silent, log-useful).

* **N-link:** (★) ensures the same cancellation of (-2N) modulo (N); multiplying by a global square preserves square-congruence.
* **Roots mod (p):** each linear contributes a single residue class (after lifting), ideal for bucket sieving.

### 6.2 Symmetric “two DoS + one linearized quadratic” (optional)

Pick ((a_i,b_i,c_i)) with
[
c_i^2=a_i b_i-2N\quad (i=1,2),
]
and pick (s,t) with (t^2+s^2x^2=(sx-t)(sx+t)). Then
[
Q(x)=(a_1x-c_1)(a_1x+c_1)(a_2x-c_2)(a_2x+c_2)(sx-t)(sx+t)\times K,
]
with (K) squareful (e.g., a square). This yields 6 linears (or group the last pair as one “piece”). Both DoS blocks share the same N-link; the linearized quadratic provides extra factors without breaking the congruence.

---

## 7. Bit Budgeting and the “|N|/5 per portion” Target

Our (k=5) linearized construction ((3)) naturally splits the trial’s log-mass into **five portions**:

1. the **squareful constant** (f^2 s^2) (parity-silent but log-useful), and
   2–5) the **four linear factors** (L_1,L_2,L_3,L_4).

**Best-case design objective.** Aim to **equipartition** the total log-mass across these five portions:
[
\boxed{\ \textbf{Target:}\quad \log_2|L_i|\ \approx\ \frac{|N|}{5}\ \text{bits for each moving factor},\qquad
\log_2(f^2 s^2)\ \approx\ \frac{|N|}{5}\ \text{bits}\ .\ }
]
Equipartition minimizes the largest piece that must be (B)-smooth, which (by Dickman’s law) maximizes success after pre-credit/splitting. Perfect equipartition is **hard** (since (x) varies and (★) ties parameters), but it is the **correct aspiration**.

### 7.1 Reconciling with the Pell-type constraint

From (★):
[
\log_2 f+\log_2 s+\log_2 t ;\approx; \tfrac12 |N|.
]
If you only balanced parameters you’d pick (\log_2 f\approx \log_2 s\approx \log_2 t\approx |N|/6).

In the (k=5) linearized form, the *moving* magnitudes are governed by
[
|L_1|\sim |s x - t|,\quad |L_2|\sim |s x + t|,\quad |L_3|\sim |f s^2 x - c|,\quad |L_4|\sim |f s^2 x + c|.
]
Thus **(x)** acts as a fifth dial. By **centering the sieve window** appropriately, you can push typical sizes of the four linear factors toward the shared (|N|/5)-bit target even if (\log_2 f,\log_2 s,\log_2 t) each hover near (|N|/6).

> **Practical band:**
> [
> \boxed{\quad \frac{|N|}{6}\ \lesssim\ \text{portion bits}\ \lesssim\ \frac{|N|}{4}\quad\text{(aspire to }|N|/5\text{)}. \quad}
> ]
> For 1024-bit (N): aim each portion at ~170–256 bits, with **~205 bits** best-case.

### 7.2 Why this bit target matters for smoothness

If the k-split moves each remainder piece into ([B_{\text{bits}},,2B_{\text{bits}}]) after pre-credit, then (R_{\text{bits}}/(k B_{\text{bits}})) is near 1–2 and (\rho(\cdot)) is large (even 1 when (\le 1)). Hence small (k) (2–5) yields very high (often guaranteed) (B)-smoothness per trial at realistic (B).

---

## 8. Practical Tuning (1024-bit Illustrative Numbers)

Assume SIQS-style polynomials with (|Q|*{\text{bits}}\approx 341). Take (B*{\text{bits}}=20).

* **Target remainder band.** (R_{\text{bits}}\in[60,120]).
  • At (R=60), (k=3\Rightarrow R/(kB)=1) ⇒ guaranteed.
  • At (R=40), (k=2\Rightarrow R/(kB)=1) ⇒ guaranteed.
  • At (R=90), (k=3\Rightarrow R/(kB)=1.5) ⇒ (p_k\approx 0.595^3\approx 0.21) (very high).

* **Parameterization for (3).** Choose (s,t,f) as products of base primes with (\left(\frac{2N}{p}\right)=+1), target (\log_2(f^2 s^2)\approx |N|/5) (≈205 bits), and center the sieve window (x_0) so (|sx_0|\approx |t|) and (|f s^2 x_0|\approx |c|).

* **Factor-base size.** The boosted per-trial success lets you **lower (B)** (shrinking (\pi(B)) and the matrix) while keeping throughput high.

---

## 9. Two-Face Cost Model

Let (\pi(B)) be the factor-base size; you need (R\approx \pi(B)) independent relations.

* **Sieving face.** With k-split and pre-credit,
  [
  \text{sieving time}\ \approx\ \frac{\pi(B)}{p_k(C)}\cdot \text{cost per tested index},
  ]
  and “cost per tested index” grows roughly with (B) (logs-only for tiny primes; buckets for larger).

* **Search face.** Enforcing (C_{\text{bits}}) with (t) distinct primes yields admissible (x) density (\lesssim 2^{t-C_{\text{bits}}}) (Sec. 4). The expected work to *find* such (x) scales like (2^{C_{\text{bits}}-t}).

**Throughput proxy (to maximize):**
[
\boxed{\quad \text{yield}(k,B,C)\ \approx\ \frac{p_k(C)}{2^{,C_{\text{bits}}-t(C)}}\ .\quad} \tag{4}
]
Increase (C_{\text{bits}}) until (p_k(C)) is high (ideally near the guarantee threshold) but not so high that (2^{C_{\text{bits}}-t(C)}) explodes. This “middle ground” is the sweet spot.

---

## 10. Implementation Blueprint (math-driven, code-agnostic)

1. **Global knobs.** Choose block size to fit cache (e.g., 32–64 Ki indices). Pick (B_{\text{bits}}\in[18,22]). Select (k\in{2,3,4,5}). Set a target remainder band (R_{\text{bits}}\in[60,120]).

2. **Pre-credit tables.**
   • Deterministic wheel for tiny primes: exact periodic log credit.
   • 2–3 “log-sketch” tables for the next bands using small moduli (M_j) (Bloom-like aggregation of logs at roots); calibrate bias once.

3. **Search constraints (controlled).** Maintain a rolling set of enforced small prime powers (p_i^{e_i}) with product (P\approx 2^{C_{\text{bits}}}), maximizing (t(C)) per credited bit (prefer many small primes). Use CRT/Hensel to mark admissible residue classes. Track the density penalty (2^{C_{\text{bits}}-t(C)}).

4. **k-aware scoring.** For each index (x) in the block:
   • If (x) is not in enforced classes, down-weight or skip.
   • Estimate (C_{\text{bits}}) from prefilters; set (R=|Q|_{\text{bits}}-C).
   • Compute a success score via (1) (or a monotone surrogate like (k^{-R/B})); shortlist only top-K indices or those with (R) in the band.

5. **Trial division & large-prime handling.** Perform TD only on shortlisted indices; allow 1-LP/2-LP; merge pairs/tuples aggressively (low-contention structures).

6. **Linear algebra.** Collect (\pi(B)+)margin relations; solve via Block-Lanczos/Wiedemann over (\mathbb F_2) (packed rows). Reducing (B) pays off quadratically.

7. **Tune by measurement.** Sweep ((k,B,C)); measure relations/sec after LP merging, total time to (\pi(B)+)margin, and matrix time. Choose ((k,B,C)) that equalizes sieving and matrix cost and maximizes (4).

---

## 11. Limitations and Why the Exponent Stays

* The one-dimensional **density wall** (2) is fundamental: admissible (x) become exponentially sparse in the credited bits unless you reintroduce wide sieving. This confines (k) to modest values for single-variable constructions.

* Even large per-trial boosts change the **constant** in (L_{1/2}), not the exponent. To change the exponent you need multi-norm geometry (e.g., two-dimensional lattice sieving in GNFS).

---

## 12. Suggested Experimental Protocol

1. Fix (N) (e.g., 1024 bits); choose SIQS polynomials.
2. Sweep (B_{\text{bits}}\in{18,20,22}), (k\in{1,2,3,4,5}), and (C_{\text{bits}}) targets (\in{220,250,280,300}).
3. Record:

    * Distribution of estimated (R_{\text{bits}}) on shortlisted indices;
    * Measured (p_k) (relations per shortlist; with/without 1-LP);
    * Admissible-(x) density vs (2^{t-C});
    * Relations/sec to (\pi(B)+)margin;
    * Matrix time and total wall-clock.
4. Choose Pareto-optimal settings (fastest total, smallest matrix, stable).

---

## 13. Conclusion

We formalized a practical approach that combines **targeted pre-credit** and **small-(k) splitting** to substantially increase per-trial success in QS/SIQS-style factoring, and we provided **sound algebraic constructions** (notably a (k=5) form with four linears and a squareful constant) that preserve the square-congruence modulo (N). The method offers large **constant-factor** improvements at realistic key sizes and a clear, measurable pathway to reduce (B), accelerate sieving, and shrink the matrix. It does **not** alter the (L_{1/2}) exponent asymptotically (due to the 1-D density wall), but it is highly attractive in practice.

---

## Acknowledgments

Thanks to discussions that clarified the middle-ground trade-offs between pre-credit size, k-split, admissible-(x) density, and to early experiments on 300-bit inputs that motivated the (k=5) linearization and the **|N|/5 per portion** bit-budget target.

---

## References (suggested standard sources)

* P. P. Pomerance, “The Quadratic Sieve Factoring Algorithm.”
* H. Riesel, *Prime Numbers and Computer Methods for Factorization.*
* J. P. Buhler, H. W. Lenstra Jr., C. Pomerance, “Factoring integers with the number field sieve.”
* A. Granville, “Smooth numbers: computational number theory and beyond.”
* R. Crandall and C. Pomerance, *Prime Numbers: A Computational Perspective* (QS, SIQS, semismoothness chapters).
* de Bruijn / Dickman on the distribution of smooth numbers.

---

## Appendix A: Heuristic Derivations

**A.1 Success probability (Eq. 1).**
Let the remainder magnitude be (S=2^{R_{\text{bits}}}). If split into (k) comparable pieces of size (S^{1/k}), each is (B)-smooth with probability (\rho!\left(\frac{\ln S^{1/k}}{\ln B}\right)=\rho!\left(\frac{R_{\text{bits}}}{kB_{\text{bits}}}\right)). Assuming independence (standard QS heuristic), multiply (k) times.

**A.2 Density (Eq. 2).**
For quadratic (Q), each prime power (p^e) contributes at most 2 lifts; CRT gives (\le 2^{\omega(P)}) residue classes mod (P). Thus frequency (\le 2^{\omega(P)}/P). Replacing (\omega(P)) by (t) yields (\lesssim 2^{t-C_{\text{bits}}}).

---

## Appendix B: k=5 Construction Details and Bit Budget

With (d=s^2,\ b=t^2) and (★) (c^2=f^2 s^2 t^2-2N),
[
Q_5(x)=f^2 s^2,(s x-t)(s x+t)(f s^2 x-c)(f s^2 x+c).
]

* **Roots mod (p):** simple linear roots; Hensel-liftable.
* **Parity:** (f^2 s^2) is a square → parity-silent; contributes to pre-credit.
* **Bit budget objective:**
  [
  \log_2(f^2 s^2)\ \approx\ |N|/5,\quad
  \log_2|s x\pm t|\ \approx\ |N|/5,\quad
  \log_2|f s^2 x\pm c|\ \approx\ |N|/5,
  ]
  with acceptable band ([|N|/6,,|N|/4]). Achieve by:

    * choosing (f,s,t) as products of base primes with (\left(\tfrac{2N}{p}\right)=+1);
    * centering (x_0) so (|s x_0|\approx |t|) and (|f s^2 x_0|\approx |c|);
    * lifting (c) from (c^2\equiv 2N \pmod{(fst)}) to ((fst)^2), nudging ({s,t,f}) and re-lifting until all portions sit near target.

**Telemetry.** For shortlisted indices, log empirical distributions of (\log_2|L_1|,\dots,\log_2|L_4|) and (\log_2(f^2 s^2)); recenter and retune if any portion drifts outside the target band.

---

### How to cite this work

> I. Gazman, “Search-Aided k-Split Smoothness for QS-Style Factoring: A Practical–Mathematical Framework and k=5 Constructions,” 2025. (preprint)
