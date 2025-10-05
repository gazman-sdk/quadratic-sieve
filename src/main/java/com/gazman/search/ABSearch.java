package com.gazman.search;

import com.gazman.factor.Logger;
import com.gazman.math.MathUtils;

import java.math.BigInteger;
import java.util.*;

/**
 * Stand-alone SIQS-style polynomial chooser:
 *  - Builds A = product of distinct primes with (N|p) = +1 near target bits (~|N|/5).
 *  - Solves B^2 ≡ N (mod A) via Tonelli + CRT (squarefree A).
 *  - Recenters B close to A*x0 for x0 in [0, xMax) and scores |Ax0±B| vs the portion band.
 *  - Emits top-K (A,B,x0) tuples with per-second speed logs.
 *
 * No dependency on factor-base or sieve; only uses MathUtils.ressol(...) for sqrt mod p.
 */
public class ABSearch extends Logger {

    public static final class Poly {
        public final BigInteger A, B;
        public final long x0;
        public final int A_bits, B_bits, Lp_bits, Lm_bits, primesCount;

        public Poly(BigInteger A, BigInteger B, long x0,
                    int A_bits, int B_bits, int Lp_bits, int Lm_bits, int primesCount) {
            this.A = A; this.B = B; this.x0 = x0;
            this.A_bits = A_bits; this.B_bits = B_bits;
            this.Lp_bits = Lp_bits; this.Lm_bits = Lm_bits; this.primesCount = primesCount;
        }
    }

    // Accept band around |N|/5
    private static final int BAND_LO = 6;
    private static final int BAND_HI = 8;

    private final BigInteger N;
    private final Random rnd = new Random(31337);

    public ABSearch(BigInteger N) { this.N = N; }

    /** Generate primes up to 'limit' with a basic sieve. */
    private static int[] primesUpTo(int limit) {
        boolean[] isComp = new boolean[limit + 1];
        ArrayList<Integer> ps = new ArrayList<>();
        for (int i = 2; i * i <= limit; i++) {
            if (!isComp[i]) for (int j = i * i; j <= limit; j += i) isComp[j] = true;
        }
        for (int i = 2; i <= limit; i++) if (!isComp[i]) ps.add(i);
        return ps.stream().mapToInt(i -> i).toArray();
    }

    /** Main entry: run for msBudget, return up to topK polynomials. */
    public List<Poly> runFixedTime(int msBudget, int topK, int primeCeil, Integer xMaxOpt) {
        final long t0 = System.currentTimeMillis();
        final long tend = t0 + msBudget;

        // 1) Build pool: primes with (N|p)=+1
        int[] primes = primesUpTo(primeCeil);
        ArrayList<Integer> pool = new ArrayList<>(8192);
        for (int p : primes) {
            if (p <= 2) continue;
            if (MathUtils.isRootInQuadraticResidues(N, BigInteger.valueOf(p))) {
                pool.add(p);
            }
        }
        if (pool.isEmpty()) {
            log("ABSearch: pool empty — no primes with (N|p)=+1 under", primeCeil);
            return List.of();
        }

        // 2) Targets
        final int nbits = N.bitLength();
        final int portion = Math.max(24, nbits / 5);
        final int A_target = portion;     // aim A near |N|/5
        final int largestPoolPrime = pool.get(pool.size() - 1);
        final int xMax = (xMaxOpt != null) ? xMaxOpt : Math.min(200_000, largestPoolPrime);
        final int xBits = Math.max(1, Integer.toBinaryString(Math.max(3, xMax)).length() - 1);

        log("ABSearch: pool", pool.size(), "| Nbits", nbits,
                "portion", portion, "| A_target", A_target,
                "| xMax", xMax, "xBits", xBits,
                "| p_max", largestPoolPrime);

        // 3) Beam search for A near target bits
        final int BEAM = 14;
        record Node(BigInteger A, ArrayList<Integer> primes, int bits) {}

        ArrayList<Node> beam = new ArrayList<>();
        beam.add(new Node(BigInteger.ONE, new ArrayList<>(), 0));

        // Keep best tuples (worst on top)
        PriorityQueue<Poly> keep = new PriorityQueue<>(Comparator
                .comparingInt((Poly p) -> bandScore(p, portion))
                .thenComparingInt(p -> Math.max(p.Lp_bits, p.Lm_bits))
                .reversed());

        // Speed logs
        long lastLog = System.currentTimeMillis();
        long attPer = 0, crtPer = 0, accPer = 0;
        int bestLpPer = Integer.MAX_VALUE, bestLmPer = Integer.MAX_VALUE;

        while (System.currentTimeMillis() < tend) {
            // Expand one round
            ArrayList<Node> next = new ArrayList<>(BEAM * 16);
            for (Node nd : beam) {
                // try a handful of random extensions
                for (int k = 0; k < 18; k++) {
                    int p = pool.get(rnd.nextInt(pool.size()));
                    if (nd.primes.contains(p)) continue;        // squarefree A
                    BigInteger A2 = nd.A.multiply(BigInteger.valueOf(p));
                    int ab = A2.bitLength();
                    if (ab > A_target + 12) continue;           // slack
                    ArrayList<Integer> list = new ArrayList<>(nd.primes.size() + 1);
                    list.addAll(nd.primes); list.add(p);
                    next.add(new Node(A2, list, ab));
                }
            }
            // also keep current nodes to be tested
            next.addAll(beam);

            // prune by closeness to A_target and number of factors
            next.sort(Comparator
                    .comparingInt((Node n) -> Math.abs(n.bits - A_target))
                    .thenComparingInt(n -> -n.primes.size()));
            beam = new ArrayList<>(next.subList(0, Math.min(BEAM, next.size())));

            // For each node, solve B^2 ≡ N (mod A), then try a few centers
            for (Node nd : beam) {
                if (nd.A.equals(BigInteger.ONE)) continue;
                attPer++;

                BigInteger Bmod = sqrtNModComposite(nd.A, nd.primes);
                if (Bmod == null) continue;
                crtPer++;

                Poly bestLocal = null;
                // Try three centers (middle, 1/3, random)
                for (int trial = 0; trial < 3; trial++) {
                    long x0 = switch (trial) {
                        case 0 -> xMax / 2L;
                        case 1 -> xMax / 3L;
                        default -> 1 + rnd.nextInt(Math.max(2, xMax - 1));
                    };
                    BigInteger target = nd.A.multiply(BigInteger.valueOf(x0));
                    BigInteger B = recenterToTarget(Bmod, nd.A, target);

                    BigInteger Lp = nd.A.multiply(BigInteger.valueOf(x0)).add(B).abs();
                    BigInteger Lm = nd.A.multiply(BigInteger.valueOf(x0)).subtract(B).abs();
                    int lpBits = Lp.bitLength();
                    int lmBits = Lm.bitLength();

                    Poly poly = new Poly(nd.A, B, x0, nd.bits, B.bitLength(),
                            lpBits, lmBits, nd.primes.size());

                    if (bestLocal == null ||
                            bandScore(poly, portion) < bandScore(bestLocal, portion) ||
                            (bandScore(poly, portion) == bandScore(bestLocal, portion)
                                    && Math.max(lpBits, lmBits) < Math.max(bestLocal.Lp_bits, bestLocal.Lm_bits))) {
                        bestLocal = poly;
                    }
                }

                if (bestLocal != null) {
                    bestLpPer = Math.min(bestLpPer, bestLocal.Lp_bits);
                    bestLmPer = Math.min(bestLmPer, bestLocal.Lm_bits);

                    if (inBand(bestLocal.Lp_bits, portion) && inBand(bestLocal.Lm_bits, portion)) {
                        accPer++;
                        if (keep.size() < topK) keep.add(bestLocal);
                        else {
                            Poly worst = keep.peek();
                            if (worseThan(worst, bestLocal, portion)) {
                                keep.poll(); keep.add(bestLocal);
                            }
                        }
                    }
                }

                // per-second speed log
                long now = System.currentTimeMillis();
                if (now - lastLog >= 1000) {
                    double sec = (now - lastLog) / 1000.0;
                    log("ABSearchSpeed:",
                            Math.round(attPer / sec) / 1000, "kAtt/s,",
                            "crt", Math.round(crtPer / sec) / 1000, "k/s,",
                            "accept", Math.round(accPer / sec), "/s |",
                            "keep", keep.size(), "| best(L+,L-)",
                            (bestLpPer == Integer.MAX_VALUE ? -1 : bestLpPer),
                            (bestLmPer == Integer.MAX_VALUE ? -1 : bestLmPer));
                    attPer = crtPer = accPer = 0;
                    bestLpPer = bestLmPer = Integer.MAX_VALUE;
                    lastLog = now;
                }

                if (System.currentTimeMillis() >= tend) break;
            }
        }

        ArrayList<Poly> out = new ArrayList<>(keep);
        out.sort(Comparator
                .comparingInt((Poly p) -> bandScore(p, portion))
                .thenComparingInt(p -> Math.max(p.Lp_bits, p.Lm_bits)));

        if (out.isEmpty()) {
            log("ABSearch: no tuples found | Nbits", N.bitLength(),
                    "portion", portion, "A_target", A_target, "xMax", xMax);
        } else {
            log("ABSearch: top", Math.min(topK, out.size()), "polynomials (closest to portion band)");
            for (int i = 0; i < Math.min(topK, out.size()); i++) {
                Poly P = out.get(i);
                log("  #", (i + 1), ": A", P.A_bits, "bits, B", P.B_bits, "bits",
                        "| x0", P.x0, "| L+~", P.Lp_bits, "L-~", P.Lm_bits,
                        "| primes=", P.primesCount);
            }
            Poly best = out.get(0);
            log("SearchSummary[SIQS]: A", best.A, "(", best.A_bits, ")",
                    "B", best.B, "(", best.B_bits, ")",
                    "x0", best.x0, "| L+~", best.Lp_bits, "L-~", best.Lm_bits,
                    "| primes=", best.primesCount);
        }
        return out;
    }

    // ---------- number helpers ----------

    private static boolean inBand(int x, int portion) {
        return x >= portion - BAND_LO && x <= portion + BAND_HI;
    }

    private static int bandScore(Poly p, int portion) {
        return Math.abs(p.Lp_bits - portion) + Math.abs(p.Lm_bits - portion);
    }

    private static boolean worseThan(Poly a, Poly b, int portion) {
        int sa = bandScore(a, portion), sb = bandScore(b, portion);
        if (sa != sb) return sa > sb;
        return Math.max(a.Lp_bits, a.Lm_bits) > Math.max(b.Lp_bits, b.Lm_bits);
    }

    /** Solve B^2 ≡ N (mod A) with A squarefree (product of distinct odd primes). */
    private static BigInteger sqrtNModComposite(BigInteger A, ArrayList<Integer> primes) {
        BigInteger r = BigInteger.ZERO, m = BigInteger.ONE;
        for (int p : primes) {
            if (p <= 2) return null;
            BigInteger P = BigInteger.valueOf(p);
            long[] roots = MathUtils.ressol(p, NmodP(P));
            long r0 = -1;
            for (long v : roots) { if (v >= 0) { r0 = v; break; } }
            if (r0 < 0) return null;
            BigInteger s = BigInteger.valueOf(r0);

            BigInteger rResidue = r.mod(m);
            BigInteger r1 = crtCombine(rResidue, m, s, P);
            BigInteger r2 = crtCombine(rResidue, m, P.subtract(s).mod(P), P);
            if (r1 == null || r2 == null) return null;

            BigInteger mNew = m.multiply(P);
            // pick representative closer to 0 to keep B magnitude controlled
            BigInteger c1 = recenterToTarget(r1, mNew, BigInteger.ZERO);
            BigInteger c2 = recenterToTarget(r2, mNew, BigInteger.ZERO);
            r = (c1.abs().compareTo(c2.abs()) <= 0) ? c1 : c2;
            m = mNew;
        }
        return m.equals(A) ? r : null;
    }

    private static long NmodP(BigInteger P) {
        return BigIntegerHolder.Ntmp.mod(P).longValue(); // see holder below
    }

    /** CRT combine (r mod m) and (x mod n), gcd=1, else null. */
    private static BigInteger crtCombine(BigInteger r, BigInteger m, BigInteger x, BigInteger n) {
        BigInteger g = m.gcd(n);
        if (!g.equals(BigInteger.ONE)) return null;
        BigInteger inv = m.modInverse(n);
        BigInteger t = x.subtract(r).mod(n);
        BigInteger coeff = t.multiply(inv).mod(n);
        return r.add(m.multiply(coeff));
    }

    /** Re-center residue r (mod M) near 'target'. */
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

    // Small static holder to let sqrtNModComposite fetch N mod p without closing over 'this'
    private static final class BigIntegerHolder {
        static BigInteger Ntmp;
    }

    /** Entrypoint-friendly helper: set static N for NmodP. */
    public static void setStaticN(BigInteger N) {
        BigIntegerHolder.Ntmp = N;
    }
}
