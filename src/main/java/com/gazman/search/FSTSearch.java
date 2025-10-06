package com.gazman.search;

import com.gazman.factor.Logger;
import com.gazman.math.SqrRoot;

import java.math.BigInteger;
import java.util.*;

/**
 * FSTSearch with shape-aware preselection:
 *  - Stage 1 uses residual lower bound to gate work,
 *    then DOES a shallow prime-only factor & greedy F/S/T to predict RMS.
 *  - Stage 2 unchanged in spirit, but with an early break in trial division.
 *  - Keeps all previous speed wins (incremental sqrt, prime cache, bounded heaps).
 */
public class FSTSearch extends Logger {

    // ---------- Tunables ----------
    public static final int SEARCH_TIME_MS = 20_000;

    // Stage 2 (finalization) full trial-division limit:
    private static final int TRIAL_DIV_LIMIT = 1_000_000;

    // Stage 1 (preselection) VERY shallow limit (cheap):
    private static final int PRE_TRIAL_DIV_LIMIT = 100_000;

    private static final int MAX_FINAL_CANDIDATES = 64;
    private static final int PRESELECT_CAP = 256;

    // Weight of N-link residual (bits) added to RMS
    private static final double RESIDUAL_WEIGHT = 0.15;

    // ---------- Inputs / derived ----------
    private final BigInteger N;
    private final BigInteger FOUR_N;
    private final int nBits;
    private final int portionBits;

    // Prime sieve cache
    private static final int[] PRIMES;             // primes up to TRIAL_DIV_LIMIT
    private static final BigInteger[] PRIMES_BI;   // cached BigInteger versions

    static {
        PRIMES = sievePrimesUpTo(TRIAL_DIV_LIMIT);
        PRIMES_BI = new BigInteger[PRIMES.length];
        for (int i = 0; i < PRIMES.length; i++) {
            PRIMES_BI[i] = BigInteger.valueOf(PRIMES[i]);
        }
    }

    public FSTSearch(BigInteger n) {
        this.N = n;
        this.FOUR_N = n.shiftLeft(2);
        this.nBits = n.bitLength();
        this.portionBits = Math.max(24, nBits / 5);
    }

    private static int[] sievePrimesUpTo(int limit) {
        boolean[] isComposite = new boolean[limit + 1];
        for (int i = 2; i * i <= limit; i++) {
            if (!isComposite[i]) {
                for (int j = i * i; j <= limit; j += i) isComposite[j] = true;
            }
        }
        int count = 0;
        for (int i = 2; i <= limit; i++) if (!isComposite[i]) count++;
        int[] primes = new int[count];
        int idx = 0;
        for (int i = 2; i <= limit; i++) if (!isComposite[i]) primes[idx++] = i;
        return primes;
    }

    private static BigInteger sqrtFloor(BigInteger n) {
        return SqrRoot.bigIntSqRootFloor(n);
    }

    // ---------- Public API ----------
    public Poly findBestPolynomial(long sieveM) {
        log("FSTSearch(k=5): start | Nbits", nBits, "| target portion", portionBits, "bits");

        final long tEnd = System.nanoTime() + SEARCH_TIME_MS * 1_000_000L;

        final BigInteger sqrtN = SqrRoot.bigIntSqRootCeil(N);
        final BigInteger m0 = sqrtN.shiftLeft(1); // ~ 2*sqrt(N)

        // Prepare incremental walks for m+delta and m-delta
        BigInteger mPlus = m0;
        BigInteger mMinus = m0;

        BigInteger c2_0 = m0.multiply(m0).subtract(FOUR_N);
        if (c2_0.signum() <= 0) c2_0 = BigInteger.ONE; // guard
        BigInteger cFloor0 = sqrtFloor(c2_0);

        BigInteger c2Plus  = c2_0;
        BigInteger cPlus   = cFloor0;

        BigInteger c2Minus = c2_0;
        BigInteger cMinus  = cFloor0;

        // Keep best PRESELECT_CAP by predicted full score (smaller is better)
        final PriorityQueue<PreCand> preHeap = new PriorityQueue<>(
                PRESELECT_CAP, Comparator.comparingDouble((PreCand p) -> p.preScore).reversed()
        );

        long checked = 0;
        long lastLog = System.nanoTime();

        // ---- Stage 1: residual-gated, shape-aware preselection ----
        for (long delta = 0; System.nanoTime() < tEnd; delta++) {

            // mPlus (delta >= 0)
            tryStashCandidate(preHeap, mPlus, cPlus, c2Plus, sieveM);

            // mMinus (delta > 0)
            if (delta > 0 && mMinus.signum() > 0) {
                tryStashCandidate(preHeap, mMinus, cMinus, c2Minus, sieveM);
            }

            checked += (delta == 0 ? 1 : 2);

            // Advance incrementally
            {   // mPlus -> mPlus+1
                BigInteger twoM = mPlus.shiftLeft(1);
                c2Plus = c2Plus.add(twoM).add(BigInteger.ONE);
                mPlus = mPlus.add(BigInteger.ONE);
                cPlus = adjustSqrtFloorUp(cPlus, c2Plus);
            }
            if (mMinus.signum() > 0) { // mMinus -> mMinus-1
                BigInteger twoM = mMinus.shiftLeft(1);
                c2Minus = c2Minus.subtract(twoM).add(BigInteger.ONE);
                mMinus = mMinus.subtract(BigInteger.ONE);
                cMinus = adjustSqrtFloorDown(cMinus, c2Minus);
            }

            long now = System.nanoTime();
            if (now - lastLog >= 1_000_000_000L) {
                lastLog = now;
                log("FSTSearch[stage1]: checked~", checked, "m's | preHeap", preHeap.size());
            }

            // Leave a couple seconds for stage 2 once we're full
            if (preHeap.size() >= PRESELECT_CAP && now + 3_000_000_000L > tEnd) break;
        }

        if (preHeap.isEmpty()) {
            log("FSTSearch: preselect empty; manufacturing fallback.");
            PreCand p = buildPreCandResidualOnly(m0, cFloor0, c2_0);
            return finalizeOne(p, sieveM);
        }

        // ---- Stage 2: finalize best few ----
        ArrayList<PreCand> bestPre = new ArrayList<>(preHeap.size());
        while (!preHeap.isEmpty()) bestPre.add(preHeap.poll());
        bestPre.sort(Comparator.comparingDouble(p -> p.preScore)); // best first

        ArrayList<Poly> finals = new ArrayList<>();
        for (int i = 0; i < bestPre.size(); i++) {
            if (System.nanoTime() >= tEnd) break;
            if (finals.size() >= MAX_FINAL_CANDIDATES) break;

            PreCand p = bestPre.get(i);
            Poly poly = finalizeOne(p, sieveM);
            finals.add(poly);

            if ((i & 7) == 7) log("FSTSearch[stage2]: finalized", finals.size(), "/", Math.min(bestPre.size(), MAX_FINAL_CANDIDATES));
        }

        if (finals.isEmpty()) {
            log("FSTSearch: finals empty; forcing finalize best pre.");
            finals.add(finalizeOne(bestPre.get(0), sieveM));
        }

        finals.sort(Comparator.comparingDouble(p -> p.score));
        Poly best = finals.get(0);
        log("FSTSearch: best score", String.format("%.2f", best.score), "| fBits", best.f.bitLength(),
                "sBits", best.s.bitLength(), "tBits", best.t.bitLength(),
                "| cBits", best.c.bitLength(), "| linkResidualBits", best.linkResidualBits);
        return best;
    }

    // ---------- Stage 1: gated preselection ----------
    private static final class PreCand {
        final BigInteger m, cFloor, c2;
        final int residualBits;
        final double preScore; // predicted full score (RMS + residual weight)

        private PreCand(BigInteger m, BigInteger cFloor, BigInteger c2, int residualBits, double preScore) {
            this.m = m;
            this.cFloor = cFloor;
            this.c2 = c2;
            this.residualBits = residualBits;
            this.preScore = preScore;
        }
    }

    private void tryStashCandidate(PriorityQueue<PreCand> heap,
                                   BigInteger m, BigInteger cFloor, BigInteger c2,
                                   long sieveM) {
        if (m.signum() <= 0 || c2.signum() <= 0) return;

        // Cheap lower bound: residual contribution only
        int residualBits = computeResidualBits(cFloor, c2);
        double lowerBound = residualBits * RESIDUAL_WEIGHT;

        // Only do shape-aware estimation if it might be competitive against worst in heap
        PreCand worst = heap.peek(); // reversed comparator -> largest at head
        boolean shouldEstimate = heap.size() < PRESELECT_CAP || (worst != null && lowerBound < worst.preScore);

        if (!shouldEstimate) return;

        // Build predicted pre-score with a very shallow factor and greedy split
        PreCand pc = buildPreCandWithShape(m, cFloor, c2, residualBits, sieveM);

        if (heap.size() < PRESELECT_CAP) {
            heap.add(pc);
        } else if (pc.preScore < heap.peek().preScore) {
            heap.poll();
            heap.add(pc);
        }
    }

    private PreCand buildPreCandResidualOnly(BigInteger m, BigInteger cFloor, BigInteger c2) {
        int residualBits = computeResidualBits(cFloor, c2);
        return new PreCand(m, cFloor, c2, residualBits, residualBits * RESIDUAL_WEIGHT);
    }

    private int computeResidualBits(BigInteger cFloor, BigInteger c2) {
        BigInteger c2Floor = cFloor.multiply(cFloor);
        if (c2.equals(c2Floor)) return 0;

        BigInteger rFloor = c2.subtract(c2Floor);
        BigInteger cCeil  = cFloor.add(BigInteger.ONE);
        BigInteger rCeil  = cCeil.multiply(cCeil).subtract(c2);
        BigInteger rBest  = (rCeil.compareTo(rFloor) < 0) ? rCeil : rFloor;
        return rBest.signum() == 0 ? 0 : rBest.bitLength();
    }

    private PreCand buildPreCandWithShape(BigInteger m, BigInteger cFloor, BigInteger c2,
                                          int residualBits, long sieveM) {
        // Choose better of floor/ceil c (no fresh sqrt)
        BigInteger c, rBest;
        BigInteger c2Floor = cFloor.multiply(cFloor);
        if (c2.equals(c2Floor)) {
            c = cFloor;
            rBest = BigInteger.ZERO;
        } else {
            BigInteger rFloor = c2.subtract(c2Floor);
            BigInteger cCeil  = cFloor.add(BigInteger.ONE);
            BigInteger rCeil  = cCeil.multiply(cCeil).subtract(c2);
            if (rCeil.compareTo(rFloor) < 0) {
                c = cCeil;
                rBest = rCeil;
            } else {
                c = cFloor;
                rBest = rFloor;
            }
        }

        // VERY shallow factoring for a shape estimate
        List<BigInteger> factors = factorByPrimeDivisionBounded(m, PRE_TRIAL_DIV_LIMIT);
        Split fst = greedySplitToFST(factors, sieveM);

        Poly approx = new Poly(fst.f, fst.s, fst.t, c, residualBits);
        approx.score(sieveM, portionBits);

        // Use the same score as Stage 2 would (RMS + residual weight)
        double preScore = approx.score;

        return new PreCand(m, cFloor, c2, residualBits, preScore);
    }

    // ---------- Stage 2: finalize (full factor, greedy, exact scoring) ----------
    private Poly finalizeOne(PreCand pc, long sieveM) {
        // choose c again (same logic)
        BigInteger cFloor = pc.cFloor;
        BigInteger c2 = pc.c2;

        BigInteger c, rBest;
        BigInteger c2Floor = cFloor.multiply(cFloor);
        if (c2.equals(c2Floor)) {
            c = cFloor;
            rBest = BigInteger.ZERO;
        } else {
            BigInteger rFloor = c2.subtract(c2Floor);
            BigInteger cCeil  = cFloor.add(BigInteger.ONE);
            BigInteger rCeil  = cCeil.multiply(cCeil).subtract(c2);
            if (rCeil.compareTo(rFloor) < 0) { c = cCeil; rBest = rCeil; }
            else { c = cFloor; rBest = rFloor; }
        }
        int residualBits = rBest.signum() == 0 ? 0 : rBest.bitLength();

        // Full factor (up to 1e6), with early break when p*p > n
        List<BigInteger> mFactors = factorByPrimeDivision(pc.m);

        Split fst = greedySplitToFST(mFactors, sieveM);
        Poly poly = new Poly(fst.f, fst.s, fst.t, c, residualBits);
        poly.score(sieveM, portionBits);
        return poly;
    }

    // ---------- Factoring ----------
    private static List<BigInteger> factorByPrimeDivisionBounded(BigInteger m, int limit) {
        ArrayList<BigInteger> res = new ArrayList<>();
        BigInteger n = m;

        // Factor 2's
        int twos = n.getLowestSetBit();
        if (twos > 0) {
            n = n.shiftRight(twos);
            for (int i = 0; i < twos; i++) res.add(BigInteger.TWO);
        }
        if (n.equals(BigInteger.ONE)) return res;

        for (int i = 1; i < PRIMES.length; i++) { // start at prime 3
            int pInt = PRIMES[i];
            if (pInt > limit) break;
            BigInteger p = PRIMES_BI[i];

            // divide out p
            while (true) {
                BigInteger[] qr = n.divideAndRemainder(p);
                if (qr[1].signum() == 0) {
                    res.add(p);
                    n = qr[0];
                    if (n.equals(BigInteger.ONE)) break;
                } else break;
            }

            // Early break if remaining n is prime (p*p > n)
            long p2 = (long) pInt * (long) pInt;
            if (BigInteger.valueOf(p2).compareTo(n) > 0) break;
        }

        if (n.compareTo(BigInteger.ONE) > 0) res.add(n);
        res.sort(Comparator.comparingInt(BigInteger::bitLength).reversed());
        return res;
    }

    private static List<BigInteger> factorByPrimeDivision(BigInteger m) {
        ArrayList<BigInteger> res = new ArrayList<>();
        BigInteger n = m;

        int twos = n.getLowestSetBit();
        if (twos > 0) {
            n = n.shiftRight(twos);
            for (int i = 0; i < twos; i++) res.add(BigInteger.TWO);
        }
        if (n.equals(BigInteger.ONE)) return res;

        for (int i = 1; i < PRIMES.length; i++) { // start at prime 3
            int pInt = PRIMES[i];
            BigInteger p = PRIMES_BI[i];

            while (true) {
                BigInteger[] qr = n.divideAndRemainder(p);
                if (qr[1].signum() == 0) {
                    res.add(p);
                    n = qr[0];
                    if (n.equals(BigInteger.ONE)) break;
                } else break;
            }

            // Early break if remaining n is prime (p*p > n)
            long p2 = (long) pInt * (long) pInt;
            if (BigInteger.valueOf(p2).compareTo(n) > 0) break;
        }

        if (n.compareTo(BigInteger.ONE) > 0) res.add(n);
        res.sort(Comparator.comparingInt(BigInteger::bitLength).reversed());
        return res;
    }

    // ---------- Greedy F/S/T ----------
    private Split greedySplitToFST(List<BigInteger> factors, long sieveM) {
        BigInteger f = BigInteger.ONE, s = BigInteger.ONE, t = BigInteger.ONE;
        if (factors.isEmpty()) return new Split(f, s, t);

        for (BigInteger g : factors) {
            double bestScore = Double.POSITIVE_INFINITY;
            int bestBin = -1;

            double scoreF = scorePreview(f.multiply(g), s, t, sieveM);
            if (scoreF < bestScore) { bestScore = scoreF; bestBin = 0; }

            double scoreS = scorePreview(f, s.multiply(g), t, sieveM);
            if (scoreS < bestScore) { bestScore = scoreS; bestBin = 1; }

            double scoreT = scorePreview(f, s, t.multiply(g), sieveM);
            if (scoreT < bestScore) { bestBin = 2; }

            if (bestBin == 0) f = f.multiply(g);
            else if (bestBin == 1) s = s.multiply(g);
            else t = t.multiply(g);
        }
        return new Split(f, s, t);
    }

    private double scorePreview(BigInteger f, BigInteger s, BigInteger t, long M) {
        Poly tmp = new Poly(f, s, t, BigInteger.ONE, 0);
        tmp.score(M, portionBits);
        return tmp.score;
    }

    private record Split(BigInteger f, BigInteger s, BigInteger t) {}

    // ---------- Incremental sqrt adjusters ----------
    private static BigInteger adjustSqrtFloorUp(BigInteger c, BigInteger c2) {
        BigInteger c1 = c.add(BigInteger.ONE);
        if (c1.multiply(c1).compareTo(c2) <= 0) {
            BigInteger c2p = c1.add(BigInteger.ONE);
            if (c2p.multiply(c2p).compareTo(c2) <= 0) return c2p;
            return c1;
        }
        if (c.multiply(c).compareTo(c2) > 0) return c.subtract(BigInteger.ONE);
        return c;
    }

    private static BigInteger adjustSqrtFloorDown(BigInteger c, BigInteger c2) {
        if (c.multiply(c).compareTo(c2) > 0) {
            BigInteger cm1 = c.subtract(BigInteger.ONE);
            if (cm1.signum() < 0) return BigInteger.ZERO;
            if (cm1.multiply(cm1).compareTo(c2) > 0) return cm1.subtract(BigInteger.ONE).max(BigInteger.ZERO);
            return cm1;
        }
        BigInteger c1 = c.add(BigInteger.ONE);
        if (c1.multiply(c1).compareTo(c2) <= 0) return c1;
        return c;
    }

    // ---------- Result type ----------
    public static class Poly {
        public final BigInteger f, s, t, c;
        public final int linkResidualBits; // 0 for exact; >0 for approximate
        public double score;

        public Poly(BigInteger f, BigInteger s, BigInteger t, BigInteger c, int linkResidualBits) {
            this.f = f;
            this.s = s;
            this.t = t;
            this.c = c;
            this.linkResidualBits = linkResidualBits;
        }

        private static double sq(double x) { return x * x; }

        /**
         * Score k=5 portions:
         * L1 = |s*M - t|,  L2 = |s*M + t|,
         * L3 = |f*s*s*M - c|,  L4 = |f*s*s*M + c|,
         * const = (f*s)^2.
         * Minimize RMS to targetBits + link residual penalty.
         */
        public void score(long M, int targetBits) {
            if (M <= 0) M = 1;
            BigInteger MB = BigInteger.valueOf(M);
            BigInteger fs2 = f.multiply(s).multiply(s);

            int l1 = s.multiply(MB).subtract(t).abs().bitLength();
            int l2 = s.multiply(MB).add(t).abs().bitLength();
            int l3 = fs2.multiply(MB).subtract(c).abs().bitLength();
            int l4 = fs2.multiply(MB).add(c).abs().bitLength();
            int cons = f.multiply(s).bitLength() * 2;

            double rms = Math.sqrt(
                    (sq(l1 - targetBits) +
                            sq(l2 - targetBits) +
                            sq(l3 - targetBits) +
                            sq(l4 - targetBits) +
                            sq(cons - targetBits)) / 5.0
            );
            this.score = rms + linkResidualBits * RESIDUAL_WEIGHT;
        }

        @Override
        public String toString() {
            return "Poly{f=" + f + ", s=" + s + ", t=" + t + ", c=" + c +
                    ", linkResidualBits=" + linkResidualBits + ", score=" + score + "}";
        }
    }
}
