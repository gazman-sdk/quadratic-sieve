# Time cut from 10 sec to 7 sec

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