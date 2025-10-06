package com.gazman.search;

import com.gazman.factor.Logger;
import com.gazman.math.MathUtils;

import java.math.BigInteger;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Stand-alone SIQS-style polynomial chooser (MULTITHREADED).
 * <p>
 * - Builds A = product of distinct primes with (N|p)=+1 near target bits (~|N|/5),
 * biased toward smaller primes to make A "richer" (more prime factors).
 * - Precomputes Tonelli roots sqrt(N) mod p for all primes in the pool.
 * - Per-thread beam exploration; for each A, creates several CRT sign-pattern variants for B,
 * tries multiple centers x0, recenters B near A*x0, and scores |Ax0±B|.
 * - Prefers tuples close to the portion band and richer A (more prime factors).
 * - Aggregates top-K across threads; prints per-second speed logs from a logger thread.
 * <p>
 * Only depends on MathUtils.ressol(...)
 */
public class ABSearch extends Logger {

    // ---------- Tunables ----------
    private static final int BAND_LO = 6;          // acceptable band around portion (|N|/5)
    private static final int BAND_HI = 8;
    private static final int DESIRED_MIN_PRIMES = 6; // prefer richer A
    private static final int RICH_WEIGHT = 1;        // penalty per missing prime
    private static final int BEAM = 18;              // per-thread beam width
    private static final int EXT_TRIES_PER_NODE = 24;
    private static final int A_SLACK_BITS = 14;      // allow A up to portion + slack
    private static final int B_VARIANTS = 6;         // CRT sign-pattern variants per A
    // bias toward smaller primes when constructing A
    private static final double SMALL_FRACTION = 0.65;   // first 65% of pool considered "small"
    private static final double SMALL_DRAW_PROB = 0.85;  // prefer drawing from the small segment
    private final BigInteger N;

    public ABSearch(BigInteger N) {
        this.N = N;
    }

    /**
     * Generate primes up to 'limit' with a basic sieve.
     */
    private static int[] primesUpTo(int limit) {
        boolean[] isComp = new boolean[Math.max(3, limit + 1)];
        ArrayList<Integer> ps = new ArrayList<>();
        for (int i = 2; i * i <= limit; i++) {
            if (!isComp[i]) for (int j = i * i; j <= limit; j += i) isComp[j] = true;
        }
        for (int i = 2; i <= limit; i++) if (!isComp[i]) ps.add(i);
        return ps.stream().mapToInt(i -> i).toArray();
    }

    private static boolean inBand(int x, int portion) {
        return x >= portion - BAND_LO && x <= portion + BAND_HI;
    }

    private static int bandScore(Poly p, int portion) {
        return Math.abs(p.Lp_bits - portion) + Math.abs(p.Lm_bits - portion);
    }

    private static int richnessPenalty(Poly p) {
        return Math.max(0, DESIRED_MIN_PRIMES - p.primesCount) * RICH_WEIGHT;
    }

    private static int compositeScore(Poly p, int portion) {
        return bandScore(p, portion) + richnessPenalty(p);
    }

    // ---------- Worker ----------

    private static boolean worseThan(Poly a, Poly b, int portion) {
        int sa = compositeScore(a, portion), sb = compositeScore(b, portion);
        if (sa != sb) return sa > sb;
        return Math.max(a.Lp_bits, a.Lm_bits) > Math.max(b.Lp_bits, b.Lm_bits);
    }

    // ---------- Helpers ----------

    private static long[] dedup(long[] a) {
        Arrays.sort(a);
        int w = 0;
        for (int i = 0; i < a.length; i++) {
            if (i == 0 || a[i] != a[i - 1]) a[w++] = a[i];
        }
        return Arrays.copyOf(a, w);
    }

    /**
     * CRT combine (r mod m) and (x mod n), gcd=1, else null.
     */
    private static BigInteger crtCombine(BigInteger r, BigInteger m, BigInteger x, BigInteger n) {
        BigInteger g = m.gcd(n);
        if (!g.equals(BigInteger.ONE)) return null;
        BigInteger inv = m.modInverse(n);
        BigInteger t = x.subtract(r).mod(n);
        BigInteger coeff = t.multiply(inv).mod(n);
        return r.add(m.multiply(coeff));
    }

    /**
     * Re-center residue r (mod M) near 'target'.
     */
    private static BigInteger recenterToTarget(BigInteger r, BigInteger M, BigInteger target) {
        BigInteger num = target.subtract(r);
        BigInteger[] dq = num.divideAndRemainder(M);
        BigInteger k = dq[0];
        BigInteger half = M.shiftRight(1);
        if (dq[1].abs().compareTo(half) > 0) {
            k = dq[1].signum() >= 0 ? k.add(BigInteger.ONE) : k.subtract(BigInteger.ONE);
        }
        return r.add(k.multiply(M));
    }

    // ----- weighted sampling over pool (prefer small primes) -----
    private static double[] buildWeights(ArrayList<Integer> pool) {
        int n = pool.size();
        int smallCut = Math.max(1, (int) Math.floor(n * SMALL_FRACTION));
        double[] w = new double[n];
        double sumSmall = 0, sumLarge = 0;

        for (int i = 0; i < smallCut; i++) {
            w[i] = 1.0 / pool.get(i);
            sumSmall += w[i];
        }
        for (int i = smallCut; i < n; i++) {
            // still allow picking large primes but with lower probability
            w[i] = (1.0 / pool.get(i)) * (1.0 - SMALL_DRAW_PROB);
            sumLarge += w[i];
        }

        double total = sumSmall + sumLarge;
        double acc = 0;
        for (int i = 0; i < n; i++) {
            acc += w[i] / total;
            w[i] = acc; // cumulative
        }
        w[n - 1] = 1.0;
        return w;
    }

    private static int drawIndex(int n, double[] cumw, ThreadLocalRandom tlr) {
        double u = tlr.nextDouble();
        int lo = 0, hi = n - 1, ans = hi;
        while (lo <= hi) {
            int mid = (lo + hi) >>> 1;
            if (u <= cumw[mid]) {
                ans = mid;
                hi = mid - 1;
            } else lo = mid + 1;
        }
        return ans;
    }

    private static void updateMin(AtomicInteger box, int val) {
        for (; ; ) {
            int cur = box.get();
            if (val >= cur) return;
            if (box.compareAndSet(cur, val)) return;
        }
    }

    /**
     * Convenience: single-thread backwards-compatible entry.
     */
    public List<Poly> runFixedTime(int msBudget, int topK, int primeCeil, Integer xMaxOpt) {
        int threads = Math.max(1, Runtime.getRuntime().availableProcessors());
        return runFixedTimeParallel(msBudget, topK, primeCeil, xMaxOpt, threads);
    }

    /**
     * Multithreaded entry: uses 'threads' workers.
     */
    public List<Poly> runFixedTimeParallel(int msBudget, int topK, int primeCeil, Integer xMaxOpt, int threads) {
        final long t0 = System.currentTimeMillis();
        final long tend = t0 + msBudget;

        // 1) Build pool: primes with (N|p) = +1
        int[] primesRaw = primesUpTo(primeCeil);
        ArrayList<Integer> pool = new ArrayList<>(8192);
        for (int p : primesRaw) {
            if (p <= 2) continue;
            if (MathUtils.isRootInQuadraticResidues(N, BigInteger.valueOf(p))) {
                pool.add(p);
            }
        }
        if (pool.isEmpty()) {
            log("ABSearch: pool empty — no primes with (N|p)=+1 under", primeCeil);
            return List.of();
        }
        final int poolMax = pool.getLast();

        // 2) Precompute sqrt(N) mod p once (Tonelli roots)
        final HashMap<Integer, PrimeRoot> rootMap = new HashMap<>(pool.size() * 2);
        for (int p : pool) {
            BigInteger P = BigInteger.valueOf(p);
            long[] rr = MathUtils.ressol(p, N.mod(P).longValue());
            long r0 = -1;
            for (long v : rr)
                if (v >= 0) {
                    r0 = v;
                    break;
                }
            if (r0 < 0) continue; // should not happen if Legendre=+1
            rootMap.put(p, new PrimeRoot(p, P, BigInteger.valueOf(r0)));
        }

        // 3) Targets and params
        final int nbits = N.bitLength();
        final int portion = Math.max(24, nbits / 5);  // ~ |N|/5
        final int xMax = (xMaxOpt != null) ? xMaxOpt : Math.min(300_000, poolMax);
        final int xBits = Math.max(1, Integer.toBinaryString(Math.max(3, xMax)).length() - 1);

        log("ABSearch: pool", pool.size(), "| Nbits", nbits,
                "portion", portion, "| A_target", portion,
                "| xMax", xMax, "xBits", xBits, "| p_max", poolMax,
                "| threads", threads);

        // 4) Weights for biased prime sampling (prefer small)
        final double[] cumw = buildWeights(pool);

        // 5) Telemetry counters (shared)
        final AtomicLong attempts = new AtomicLong(0);
        final AtomicLong crts = new AtomicLong(0);
        final AtomicLong accepts = new AtomicLong(0);
        final AtomicInteger bestLp = new AtomicInteger(Integer.MAX_VALUE);
        final AtomicInteger bestLm = new AtomicInteger(Integer.MAX_VALUE);
        final AtomicBoolean stop = new AtomicBoolean(false);

        // 6) Workers
        ExecutorService exec = Executors.newFixedThreadPool(threads);
        ArrayList<Future<PriorityQueue<Poly>>> futures = new ArrayList<>(threads);

        for (int tid = 0; tid < threads; tid++) {
            final int threadId = tid;
            futures.add(exec.submit(() ->
                    workerSearch(threadId, tend, pool, cumw, rootMap, portion, portion, xMax,
                            attempts, crts, accepts, bestLp, bestLm, topK)
            ));
        }

        // 7) Logger thread
        Thread logger = new Thread(() -> {
            long last = System.currentTimeMillis();
            long a0 = 0, c0 = 0, ac0 = 0;
            while (!stop.get()) {
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException ignored) {
                }
                long now = System.currentTimeMillis();
                long a1 = attempts.get(), c1 = crts.get(), ac1 = accepts.get();
                double sec = Math.max(1, (now - last)) / 1000.0;
                log("ABSearchSpeed:",
                        Math.round((a1 - a0) / sec) / 1000, "kAtt/s,",
                        "crt", Math.round((c1 - c0) / sec) / 1000, "k/s,",
                        "accept", (int) Math.round((ac1 - ac0) / sec), "/s |",
                        "keep", "-", "| best(L+,L-)",
                        (bestLp.get() == Integer.MAX_VALUE ? -1 : bestLp.get()),
                        (bestLm.get() == Integer.MAX_VALUE ? -1 : bestLm.get()));
                last = now;
                a0 = a1;
                c0 = c1;
                ac0 = ac1;
                if (System.currentTimeMillis() >= tend) break;
            }
        }, "ABSearch-logger");
        logger.setDaemon(true);
        logger.start();

        // 8) Collect results
        ArrayList<Poly> collected = new ArrayList<>();
        try {
            for (Future<PriorityQueue<Poly>> f : futures) {
                PriorityQueue<Poly> pq = f.get();
                collected.addAll(pq);
            }
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            stop.set(true);
            exec.shutdownNow();
            try {
                logger.join(200);
            } catch (InterruptedException ignored) {
            }
        }

        // 9) Rank and print results
        collected.sort(Comparator
                .comparingInt((Poly p) -> compositeScore(p, portion))
                .thenComparingInt(p -> Math.max(p.Lp_bits, p.Lm_bits)));

        if (collected.isEmpty()) {
            log("ABSearch: no tuples found | Nbits", N.bitLength(),
                    "portion", portion);
            return List.of();
        }

        int show = Math.min(topK, collected.size());
        log("ABSearch: top", show, "polynomials (closest to portion band; richer A preferred)");
        for (int i = 0; i < show; i++) {
            Poly P = collected.get(i);
            log("  #", (i + 1), ": A", P.A_bits, "bits, B", P.B_bits, "bits",
                    "| x0", P.x0, "| L+~", P.Lp_bits, "L-~", P.Lm_bits,
                    "| primes=", P.primesCount);
        }
        Poly best = collected.getFirst();
        log("SearchSummary[SIQS]: A", best.A, "(", best.A_bits, ")",
                "B", best.B, "(", best.B_bits, ")",
                "x0", best.x0, "| L+~", best.Lp_bits, "L-~", best.Lm_bits,
                "| primes=", best.primesCount);

        return collected.subList(0, show);
    }

    private PriorityQueue<Poly> workerSearch(
            int threadId,
            long tend,
            ArrayList<Integer> pool,
            double[] cumw,
            HashMap<Integer, PrimeRoot> rootMap,
            int portion,
            int A_target,
            int xMax,
            AtomicLong attempts,
            AtomicLong crts,
            AtomicLong accepts,
            AtomicInteger bestLp,
            AtomicInteger bestLm,
            int topK
    ) {
        final ThreadLocalRandom tlr = ThreadLocalRandom.current();

        // per-thread keep (worst on top)
        PriorityQueue<Poly> keep = new PriorityQueue<>(Comparator
                .comparingInt((Poly p) -> compositeScore(p, portion))
                .thenComparingInt(p -> Math.max(p.Lp_bits, p.Lm_bits))
                .reversed());

        ArrayList<Node> beam = new ArrayList<>();
        beam.add(new Node(BigInteger.ONE, new ArrayList<>(), 0));

        while (System.currentTimeMillis() < tend) {
            // expand beam one round
            ArrayList<Node> next = new ArrayList<>(BEAM * EXT_TRIES_PER_NODE);
            for (Node nd : beam) {
                for (int k = 0; k < EXT_TRIES_PER_NODE; k++) {
                    int pickIdx = drawIndex(pool.size(), cumw, tlr);
                    int p = pool.get(pickIdx);
                    if (nd.primes.contains(p)) continue;            // squarefree A
                    BigInteger A2 = nd.A.multiply(BigInteger.valueOf(p));
                    int ab = A2.bitLength();
                    if (ab > A_target + A_SLACK_BITS) continue;     // slack above target
                    ArrayList<Integer> list = new ArrayList<>(nd.primes.size() + 1);
                    list.addAll(nd.primes);
                    list.add(p);
                    next.add(new Node(A2, list, ab));
                }
            }
            next.addAll(beam);

            next.sort(Comparator
                    .comparingInt((Node n) -> Math.abs(n.bits - A_target))
                    .thenComparingInt(n -> -n.primes.size()));
            beam = new ArrayList<>(next.subList(0, Math.min(BEAM, next.size())));

            // For each node, build B variants and test several centers
            for (Node nd : beam) {
                if (nd.A.equals(BigInteger.ONE)) continue;
                attempts.incrementAndGet();

                // Build PrimeRoot list from cache
                ArrayList<PrimeRoot> roots = new ArrayList<>(nd.primes.size());
                boolean ok = true;
                for (int p : nd.primes) {
                    PrimeRoot pr = rootMap.get(p);
                    if (pr == null) {
                        ok = false;
                        break;
                    }
                    roots.add(pr);
                }
                if (!ok) continue;

                List<BigInteger> bResidues = buildBVariants(roots, nd.A);
                if (bResidues.isEmpty()) continue;
                crts.incrementAndGet();

                long[] centers = proposeCenters(xMax, tlr);
                Poly bestLocal = null;
                for (BigInteger Bmod : bResidues) {
                    for (long x0 : centers) {
                        BigInteger target = nd.A.multiply(BigInteger.valueOf(x0));
                        BigInteger B = recenterToTarget(Bmod, nd.A, target);
                        Poly cand = scorePoly(nd, B, x0);

                        // update best L±
                        updateMin(bestLp, cand.Lp_bits);
                        updateMin(bestLm, cand.Lm_bits);

                        bestLocal = better(bestLocal, cand, portion);
                    }
                }

                if (bestLocal != null && inBand(bestLocal.Lp_bits, portion) && inBand(bestLocal.Lm_bits, portion)) {
                    accepts.incrementAndGet();
                    if (keep.size() < topK) keep.add(bestLocal);
                    else {
                        Poly worst = keep.peek();
                        if (worseThan(worst, bestLocal, portion)) {
                            keep.poll();
                            keep.add(bestLocal);
                        }
                    }
                }

                if (System.currentTimeMillis() >= tend) break;
            }
        }

        return keep;
    }

    private Poly better(Poly cur, Poly cand, int portion) {
        if (cur == null) return cand;
        return worseThan(cur, cand, portion) ? cand : cur;
    }

    private Poly scorePoly(Node nd, BigInteger B, long x0) {
        BigInteger Ax = nd.A.multiply(BigInteger.valueOf(x0));
        BigInteger Lp = Ax.add(B).abs();
        BigInteger Lm = Ax.subtract(B).abs();
        return new Poly(nd.A, B, x0, nd.bits, B.bitLength(), Lp.bitLength(), Lm.bitLength(), nd.primes.size());
    }

    /**
     * Propose centers per A (uses per-thread RNG).
     */
    private long[] proposeCenters(int xMax, ThreadLocalRandom tlr) {
        if (xMax <= 3) return new long[]{1};
        long mid = xMax / 2L;
        long oneThird = xMax / 3L;
        long twoThird = (2L * xMax) / 3L;
        long oneQuarter = xMax / 4L;
        long threeQuarter = (3L * xMax) / 4L;
        long rand = 1 + tlr.nextInt(Math.max(2, xMax - 1));
        return dedup(new long[]{mid, oneThird, twoThird, oneQuarter, threeQuarter, 1, rand});
    }

    // ----- Build several residues B (mod A) by flipping per-prime signs -----
    private List<BigInteger> buildBVariants(List<PrimeRoot> roots, BigInteger A) {
        ArrayList<BigInteger> out = new ArrayList<>(ABSearch.B_VARIANTS);
        HashSet<BigInteger> seen = new HashSet<>(ABSearch.B_VARIANTS * 2);

        // base (+ + + ...)
        BigInteger base = crtCombineAll(roots, null);
        if (base != null) {
            BigInteger canon = base.mod(A);
            seen.add(canon);
            out.add(canon);
            BigInteger neg = A.subtract(canon).mod(A);
            if (seen.add(neg)) out.add(neg);
        }

        // random patterns
        int n = roots.size();
        int flips = Math.max(1, (int) Math.ceil(n / 3.0)); // flip ~1/3 per variant
        ThreadLocalRandom tlr = ThreadLocalRandom.current();
        for (int t = 0; t < ABSearch.B_VARIANTS && out.size() < ABSearch.B_VARIANTS; t++) {
            boolean[] sign = new boolean[n];
            for (int i = 0; i < flips; i++) sign[tlr.nextInt(n)] = true; // true => use (P - r)
            BigInteger v = crtCombineAll(roots, sign);
            if (v == null) continue;
            BigInteger canon = v.mod(A);
            if (seen.add(canon)) out.add(canon);
        }
        return out;
    }

    /**
     * CRT combine all residues with optional sign flips. sign[i]==true uses (P_i - r_i).
     */
    private BigInteger crtCombineAll(List<PrimeRoot> roots, boolean[] sign) {
        BigInteger R = BigInteger.ZERO;
        BigInteger M = BigInteger.ONE;
        for (int i = 0; i < roots.size(); i++) {
            PrimeRoot pr = roots.get(i);
            BigInteger residue = (sign != null && sign[i]) ? pr.P.subtract(pr.r).mod(pr.P) : pr.r;
            BigInteger nxt = crtCombine(R.mod(M), M, residue, pr.P);
            if (nxt == null) return null;
            R = nxt;
            M = M.multiply(pr.P);
        }
        return R;
    }

    // ---------- Public result tuple ----------
    public record Poly(BigInteger A, BigInteger B, long x0, int A_bits, int B_bits, int Lp_bits, int Lm_bits,
                       int primesCount) {
    }

    /**
     * @param primes list of primes in A
     */ // Beam node (accessible so helpers can use it)
    record Node(BigInteger A, ArrayList<Integer> primes, int bits) {
    }

    /**
     * @param r one root (the other is P - r)
     */ // Root for a single prime
    private record PrimeRoot(int p, BigInteger P, BigInteger r) {
    }
}
