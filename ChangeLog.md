# Time cut from 10.5 sec to 7.9 sec

* **Single-pass sieve.**
  I merged your “add logs” and “build vectors” passes into one pass per block. We only scan primes once, then decide which indices are worth trial-dividing.

* **Bucket sieve (hit lists).**
  While sieving, for each index I also record *which* primes hit it. When a position looks promising, I trial-divide **only by those primes**, not by the whole factor base. This slashes trial-division work.

* **Thread-local batching.**
  Each worker collects B-smooths and single-large-prime relations in thread-local lists and merges once per block. No synchronization in the hot loop and far fewer allocations.

* **Fixed-point (integer) logs.**
  Precompute `log(p)` scaled by 256 and accumulate in an `int[]`. This keeps the sieve inner loop in integer math and avoids calling `Math.log` per hit.

* **Exact candidate check, then TD.**
  For each index, compute the expected `ln|Q(t)| = ln(t*(t+2*root))` once. If `expected – accumulated ≤ ln(B_max²)` (plus a tiny epsilon), *then* trial-divide. This recreates your original acceptance criterion, but only for shortlisted indices.

* **Large-prime policy = your original.**
  After trial division, accept a single large prime iff the remainder `≤ (B_max)²` and is probable prime. This restored the “Big primes found” stream that your shrinker expects.

# What changed (code structure)

* **New thread-local workspace:**
  `acc[]` (scaled log sums), `head[]/who[]/next[]` (per-index buckets as a compact linked list), and small per-thread vectors for results.

* **Sieve loop:**
  `w.prepareToMove(); while (w.testMove()) { idx = w.move(); acc[idx]+=logP; push prime i into bucket[idx]; }`

* **Candidate scan:**
  For each index, compute `expectedLn`, compare to `acc[idx]`. If inside budget, call `trialDivideBucketed(t, head[idx], …)` which:

    * builds the parity BitSet by dividing only by the primes in the bucket list (including powers),
    * emits B-smooth if remainder == 1,
    * emits single-large-prime relation if remainder ≤ (B_max)² and prime.

* **Merging + logging:**
  Results are merged once per block; the matrix/solution path is unchanged.

# Why wall-time improved but “kValues/s” dropped

Your “kValues/s” counter increments **once per block**, not per candidate processed. With the bucket sieve we spend a bit more time per block doing smarter candidate selection and far less pointless trial division, so:

* blocks/sec ↓ (your metric),
* **relations/sec ↑**, and
* total wall-time ↓ (the number you care about).

If you want that metric to reflect real throughput, add a “relations/sec” stat (B-smooths + usable large-prime relations).

# Safe tunables I left on the table (easy extra speed)

* **Two-stage early abort:** after dividing the first few bucket primes, compare a conservative log of the remainder with the remaining budget; bail early if it can’t fit. This typically trims another 20–40% of TD time.
* **Bigger block (segment) size:** making `sieveVectorBound` a larger fixed segment (e.g., 128k–512k) increases cache-friendly linear passes and reduces per-block overhead.
* **Pre-sieving small prime powers (2, 3, 5):** subtract their logs in tight unrolled loops before the general sieve.


# Time cut from 7.9 sec to 6.1 sec

awesome—glad it moved the needle. here’s a tight summary of what we changed and why it helped:

# What changed

1. Two-Large-Prime (2LP) relations with instant merging

* Added `BigPrimePairs` to collect partial relations keyed by an unordered pair `(p, q)`.
* When a second relation with the same `(p, q)` arrives, we XOR their parity vectors and multiply their `(x,y)` to produce a **B-smooth** relation immediately.
* Result: many near-smooths convert to full relations sooner → fewer candidates needed to hit the matrix target → less sieving time.

2. Killed per-index `Math.log` in the hot loop (block-constant threshold)

* Precompute a **single** log threshold per block: `LN_2ROOT_SCALED + log(t0)` instead of `log(Q(t))` for every index.
* Kept the bucketed trial division, so occasional false positives are filtered cheaply.
* Result: removes millions of transcendental calls per run → big pure CPU win; your speed jumped to ~200k values/s.

3. Kept the fast, cache-friendly sieve structure

* Still bucket by index, still only trial-divide by primes that actually hit that index (from the bucket list).
* Threads batch their finds and merge in bulk (no extra contention).

# Code touchpoints (Java)

* New class: `BigPrimePairs` (thread-safe hashmap keyed by `"min#max"`), merges equal 2LP pairs into B-smooth relations.
* `QuadraticThieve`:

    * Added fields: `LN_2ROOT_SCALED` and `BigPrimePairs bigPrimePairs`.
    * Constructor: initialize both.
    * Replaced `sieveOnePass(...)`: uses **block-constant** expected log; no per-index `Math.log`.
    * Replaced `trialDivideBucketed(...)`:

        * unchanged bucketed TD for factor-base primes,
        * single-LP path kept,
        * **new 2LP path**: cheap 64-bit split for `rem ≤ B_max^2`, verify both factors prime, then try to merge via `bigPrimePairs`.
    * Added `factorSemiprimeLE1e10(long n)` helper (trial divides up to 100k; safe because your large-prime cutoff is `B_max^2`).

# Why it sped things up

* 2LP merging increases the **yield per candidate**, so you reach ~B relations sooner.
* Removing per-index logs reduces hot-loop cost substantially; CPU now spends cycles on real sieving/TD, not transcendental math.
* Matrix phase unchanged (already ~1s), but it kicks in earlier thanks to more/faster relations.

# Sanity/guardrails

* Thread-safety: `BigPrimePairs.add(...)` is `synchronized`; per-thread batches still used for low contention elsewhere.
* Correctness: we only create a B-smooth from 2LP when we see **the same (p, q)** twice; parity vectors XORed; `x,y` multiplied.
* Bounds: the 64-bit splitter is only called when `rem ≤ B_max^2`; at least one factor ≤ 1e5, so it’s fast.