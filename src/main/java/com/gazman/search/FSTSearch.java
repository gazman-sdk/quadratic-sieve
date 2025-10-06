package com.gazman.search;

import com.gazman.factor.Logger;
import com.gazman.math.SqrRoot;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Search for k=5 coefficients (f, s, t, c) satisfying (preferably)
 * c^2 = (f s t)^2 - 4N.
 * <p>
 * Notes:
 * - Exact hits exist for semiprimes only at m = p + q (Fermat); they can be very far
 * from 2*sqrt(N), so a short linear scan usually won't find them.
 * - This implementation still checks for exact squares, but ALSO accepts approximate
 * candidates by rounding c ≈ sqrt(m^2 - 4N) and adding a small "link residual"
 * penalty to the score: residual = |c^2 - (m^2 - 4N)|.
 * - Result: you always get the best candidate within the time budget (never throw).
 */
public class FSTSearch extends Logger {

    // Tunables
    public static final int SEARCH_TIME_MS = 20_000;
    private static final int TRIAL_DIV_LIMIT = 1_000_000;
    private static final int MAX_CANDIDATES = 64;
    // Weight for the N-link residual in the score: add linkResidualBits * RESIDUAL_WEIGHT to RMS size error.
    private static final double RESIDUAL_WEIGHT = 0.15;
    private final BigInteger N;
    private final BigInteger FOUR_N;
    private final int nBits;
    private final int portionBits;

    public FSTSearch(BigInteger n) {
        this.N = n;
        this.FOUR_N = n.shiftLeft(2);
        this.nBits = n.bitLength();
        this.portionBits = Math.max(24, nBits / 5);
    }

    private static BigInteger sqrtIfPerfect(BigInteger n) {
        BigInteger r = SqrRoot.bigIntSqRootFloor(n);
        return r.multiply(r).equals(n) ? r : null;
    }

    private static List<BigInteger> factorByTrialDivision(BigInteger m, int limit) {
        ArrayList<BigInteger> res = new ArrayList<>();
        BigInteger n = m;

        // Factor 2's
        int twos = n.getLowestSetBit();
        if (twos > 0) {
            n = n.shiftRight(twos);
            for (int i = 0; i < twos; i++) res.add(BigInteger.TWO);
        }

        // Odd trial division
        for (long d = 3; d <= limit; d += 2) {
            BigInteger D = BigInteger.valueOf(d);
            BigInteger D2 = D.multiply(D);
            if (D2.compareTo(n) > 0) break;
            while (n.mod(D).signum() == 0) {
                res.add(D);
                n = n.divide(D);
            }
        }

        if (n.compareTo(BigInteger.ONE) > 0) res.add(n);
        // Sort largest factors first to help greedy balancing
        res.sort(Comparator.comparingInt(BigInteger::bitLength).reversed());
        return res;
    }

    public Poly findBestPolynomial(long sieveM) {
        log("FSTSearch(k=5): start | Nbits", nBits, "| target portion", portionBits, "bits");

        final long tend = System.currentTimeMillis() + SEARCH_TIME_MS;
        final BigInteger sqrtN = SqrRoot.bigIntSqRootCeil(N);
        final BigInteger m0 = sqrtN.shiftLeft(1); // ~ 2*sqrt(N)

        final ArrayList<Poly> bag = new ArrayList<>();
        long checked = 0;
        long lastLog = System.currentTimeMillis();

        for (long delta = 0; System.currentTimeMillis() < tend; delta++) {
            // m = m0 + delta
            if (System.currentTimeMillis() >= tend) break;
            BigInteger mPlus = m0.add(BigInteger.valueOf(delta));
            maybeRecordCandidate(mPlus, bag, sieveM);
            checked++;

            // m = m0 - delta (if positive)
            if (delta > 0) {
                if (System.currentTimeMillis() >= tend) break;
                BigInteger mMinus = m0.subtract(BigInteger.valueOf(delta));
                if (mMinus.signum() > 0) {
                    maybeRecordCandidate(mMinus, bag, sieveM);
                    checked++;
                }
            }

            long now = System.currentTimeMillis();
            if (now - lastLog >= 1000) {
                lastLog = now;
                log("FSTSearch: checked~", checked, "m's | kept", bag.size(), "candidates");
            }
        }

        if (bag.isEmpty()) {
            // This should basically never happen now (m0 gives a valid approx).
            // As a super-safe fallback, manufacture a minimal candidate from m0.
            log("FSTSearch: no candidates collected; manufacturing a fallback.");
            BigInteger m = m0.max(BigInteger.ONE);
            Poly fallback = buildApproxCandidate(m, sieveM);
            return fallback;
        }

        bag.sort(Comparator.comparingDouble(p -> p.score));
        Poly best = bag.get(0);
        log("FSTSearch: best score", String.format("%.2f", best.score), "| fBits", best.f.bitLength(),
                "sBits", best.s.bitLength(), "tBits", best.t.bitLength(),
                "| cBits", best.c.bitLength(), "| linkResidualBits", best.linkResidualBits);
        return best;
    }

    private void maybeRecordCandidate(BigInteger m, ArrayList<Poly> bag, long sieveM) {
        if (m.signum() <= 0) return;

        // Prefer exact square if present; otherwise accept best rounded c with residual penalty.
        BigInteger c2 = m.multiply(m).subtract(FOUR_N);
        if (c2.signum() <= 0) return;

        // Try exact first
        BigInteger cExact = sqrtIfPerfect(c2);
        if (cExact != null) {
            Poly exact = buildCandidateWithC(m, cExact, 0, sieveM);
            bag.add(exact);
        } else {
            // Approximate: choose the better of floor/ceil roots and penalize by residual bits.
            Poly approx = buildApproxCandidate(m, sieveM);
            bag.add(approx);
        }

        if (bag.size() > MAX_CANDIDATES * 2) {
            bag.sort(Comparator.comparingDouble(p -> p.score));
            while (bag.size() > MAX_CANDIDATES) bag.remove(bag.size() - 1);
        }
    }

    private Poly buildApproxCandidate(BigInteger m, long sieveM) {
        BigInteger c2 = m.multiply(m).subtract(FOUR_N);
        BigInteger cFloor = SqrRoot.bigIntSqRootFloor(c2);
        BigInteger rFloor = c2.subtract(cFloor.multiply(cFloor)); // >= 0
        BigInteger cBest = cFloor;
        BigInteger rBest = rFloor;

        if (!rFloor.equals(BigInteger.ZERO)) {
            BigInteger cCeil = cFloor.add(BigInteger.ONE);
            BigInteger rCeil = cCeil.multiply(cCeil).subtract(c2); // >= 0
            if (rCeil.compareTo(rFloor) < 0) {
                cBest = cCeil;
                rBest = rCeil;
            }
        }

        int residualBits = rBest.equals(BigInteger.ZERO) ? 0 : rBest.bitLength();
        return buildCandidateWithC(m, cBest, residualBits, sieveM);
    }

    private Poly buildCandidateWithC(BigInteger m, BigInteger c, int linkResidualBits, long sieveM) {
        List<BigInteger> mFactors = factorByTrialDivision(m, TRIAL_DIV_LIMIT);
        Split fst = greedySplitToFST(mFactors, sieveM);
        Poly poly = new Poly(fst.f, fst.s, fst.t, c, linkResidualBits);
        poly.score(sieveM, portionBits);
        return poly;
    }

    private Split greedySplitToFST(List<BigInteger> factors, long sieveM) {
        BigInteger f = BigInteger.ONE, s = BigInteger.ONE, t = BigInteger.ONE;
        if (factors.isEmpty()) return new Split(f, s, t);

        for (BigInteger g : factors) {
            double bestScore = Double.POSITIVE_INFINITY;
            int bestBin = -1;

            double scoreF = scorePreview(f.multiply(g), s, t, sieveM);
            if (scoreF < bestScore) {
                bestScore = scoreF;
                bestBin = 0;
            }

            double scoreS = scorePreview(f, s.multiply(g), t, sieveM);
            if (scoreS < bestScore) {
                bestScore = scoreS;
                bestBin = 1;
            }

            double scoreT = scorePreview(f, s, t.multiply(g), sieveM);
            if (scoreT < bestScore) {
                bestBin = 2;
            }

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

    private record Split(BigInteger f, BigInteger s, BigInteger t) {
    }

    // ----- Result -----
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

        private static double sq(double x) {
            return x * x;
        }

        /**
         * Score the five portions for k=5:
         * L1 = |s*M - t|,  L2 = |s*M + t|,
         * L3 = |f*s*s*M - c|,  L4 = |f*s*s*M + c|,
         * const = (f*s)^2.
         * We minimize RMS error vs 'targetBits', plus a small penalty for linkResidualBits.
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

            // Penalize deviation from exact N-link; exact hits (0 bits) win.
            this.score = rms + linkResidualBits * RESIDUAL_WEIGHT;
        }

        @Override
        public String toString() {
            return "Poly{f=" + f + ", s=" + s + ", t=" + t + ", c=" + c +
                    ", linkResidualBits=" + linkResidualBits + ", score=" + score + "}";
        }
    }
}
