package com.gazman.search;

import com.gazman.factor.Logger;
import com.gazman.math.SqrRoot;

import java.math.BigInteger;
import java.util.*;

/**
 * FSTSearch with:
 *  - QS-style residue wheel presieve for m (dramatically fewer candidates)
 *  - Shape-aware preselection using an edge-aware "elite" score
 *  - Final logs: sBits, fsBits, slopes, and line bit-sizes at x=1, 2^10, 2^20, 2^30 (both +x and -x)
 *  - Predicted edge-bit summaries for B = 2^20 and B = 2^30
 *
 * Public API preserved: findBestPolynomial(long sieveM)
 */
public class FSTSearch extends Logger {

    // ---------- Tunables ----------
    public static final int SEARCH_TIME_MS = 20_000;

    // Trial division limits (bounded prime-only)
    private static final int TRIAL_DIV_LIMIT = 1_000_000;     // Stage-2 (finalization)
    private static final int PRE_TRIAL_DIV_LIMIT = 100_000;   // Stage-1 (preselection)

    // Candidate counts
    private static final int MAX_FINAL_CANDIDATES = 64;
    private static final int PRESELECT_CAP = 256;

    // QS wheel primes (product W ≈ 255,255) — prune ~1/2 per prime.
    private static final int[] WHEEL_PRIMES = {3, 5, 7, 11, 13, 17};

    // Stage-1 edge window assumption (half width) for pre-score
    private static final int PRESELECT_B_BITS = 20;

    // Weight of the N-link residual penalty (in bits) added to RMS
    private static final double RESIDUAL_WEIGHT = 0.15;

    // ---------- Inputs / derived ----------
    private final BigInteger N;
    private final BigInteger FOUR_N;
    private final int nBits;
    private final int portionBits;

    // Prime cache up to TRIAL_DIV_LIMIT (for factoring)
    private static final int[] PRIMES;
    private static final BigInteger[] PRIMES_BI;

    static {
        PRIMES = sievePrimesUpTo(TRIAL_DIV_LIMIT);
        PRIMES_BI = new BigInteger[PRIMES.length];
        for (int i = 0; i < PRIMES.length; i++) PRIMES_BI[i] = BigInteger.valueOf(PRIMES[i]);
    }

    public FSTSearch(BigInteger n) {
        this.N = n;
        this.FOUR_N = n.shiftLeft(2);
        this.nBits = n.bitLength();
        this.portionBits = Math.max(24, nBits / 5);
    }

    // ---------- Public API ----------
    public Poly findBestPolynomial(long sieveM) {
        log("FSTSearch(k=5): start | Nbits", nBits, "| target portion", portionBits, "bits");

        final long tEnd = System.nanoTime() + SEARCH_TIME_MS * 1_000_000L;

        // Build QS residue wheel (depends on N)
        final Wheel wheel = Wheel.build(N, WHEEL_PRIMES);
        log("Wheel: W=", wheel.W, "| allowed residues=", wheel.allowedIdx.length,
                "| avg skip≈", String.format("%.1f", (double) wheel.W / Math.max(1, wheel.allowedIdx.length)));

        // m0 ≈ 2*sqrt(N)
        final BigInteger sqrtN = SqrRoot.bigIntSqRootCeil(N);
        final BigInteger m0 = sqrtN.shiftLeft(1);

        // Seed scan states on the wheel
        ScanState plus  = wheel.seedUpFrom(m0);
        ScanState minus = wheel.seedDownFrom(m0);

        // Stage-1: keep best PRESELECT_CAP by predicted score
        final PriorityQueue<PreCand> preHeap = new PriorityQueue<>(
                PRESELECT_CAP, Comparator.comparingDouble((PreCand p) -> p.preScore).reversed()
        );

        long sampled = 0;
        long lastLog = System.nanoTime();

        // ---- Stage 1: presieved scan with edge-aware preselection ----
        while (System.nanoTime() < tEnd) {
            // step up
            if (plus != null) {
                handleCandidate(preHeap, plus.m, sieveM);
                sampled++;
                plus = wheel.stepUp(plus);
            }
            if (System.nanoTime() >= tEnd) break;
            // step down
            if (minus != null) {
                handleCandidate(preHeap, minus.m, sieveM);
                sampled++;
                minus = wheel.stepDown(minus);
            }

            long now = System.nanoTime();
            if (now - lastLog >= 1_000_000_000L) {
                lastLog = now;
                log("FSTSearch[stage1]: sampled~", sampled, "m's | preHeap", preHeap.size());
            }
            // leave a few seconds for stage 2 once heap is full
            if (preHeap.size() >= PRESELECT_CAP && now + 3_000_000_000L > tEnd) break;

            // If one direction exhausted, keep going in the other
            if (plus == null && minus == null) break;
        }

        if (preHeap.isEmpty()) {
            // Fallback: evaluate m0 directly
            BigInteger c2 = m0.multiply(m0).subtract(FOUR_N);
            BigInteger cFloor = c2.signum() > 0 ? SqrRoot.bigIntSqRootFloor(c2) : BigInteger.ONE;
            PreCand pc = buildPreCandResidualOnly(m0, cFloor, c2);
            return finalizeAndLog(pc, sieveM);
        }

        // ---- Stage 2: finalize best few ----
        ArrayList<PreCand> bestPre = new ArrayList<>(preHeap.size());
        while (!preHeap.isEmpty()) bestPre.add(preHeap.poll());
        bestPre.sort(Comparator.comparingDouble(p -> p.preScore));

        ArrayList<Poly> finals = new ArrayList<>();
        for (int i = 0; i < bestPre.size() && System.nanoTime() < tEnd && finals.size() < MAX_FINAL_CANDIDATES; i++) {
            if ((i & 7) == 0) log("FSTSearch[stage2]: finalized", finals.size(), "/", Math.min(bestPre.size(), MAX_FINAL_CANDIDATES));
            finals.add(finalizeOne(bestPre.get(i), sieveM));
        }
        if (finals.isEmpty()) finals.add(finalizeOne(bestPre.get(0), sieveM));

        finals.sort(Comparator.comparingDouble(p -> p.score));
        Poly best = finals.get(0);

        // Summary line
        log("FSTSearch: best score", String.format("%.2f", best.score), "| fBits", best.f.bitLength(),
                "sBits", best.s.bitLength(), "tBits", best.t.bitLength(),
                "| cBits", best.c.bitLength(), "| linkResidualBits", best.linkResidualBits);

        // Detailed behavior logs (center, x=1, x=2^10, x=2^20, x=2^30; both +x and -x)
        logEliteSummary(best, sieveM);

        // Predicted edge bits using slope model for typical windows
        logPredictedEdgeBits(best, sieveM, 20);
        logPredictedEdgeBits(best, sieveM, 30);

        return best;
    }

    // ---------- Stage 1 helpers ----------
    private void handleCandidate(PriorityQueue<PreCand> heap, BigInteger m, long sieveM) {
        BigInteger c2 = m.multiply(m).subtract(FOUR_N);
        if (c2.signum() <= 0) return;

        BigInteger cFloor = SqrRoot.bigIntSqRootFloor(c2);
        int residualBits = computeResidualBits(cFloor, c2);
        double lowerBound = residualBits * RESIDUAL_WEIGHT;

        PreCand worst = heap.peek();
        boolean worthEstimating = heap.size() < PRESELECT_CAP || (worst != null && lowerBound < worst.preScore);
        if (!worthEstimating) return;

        // Shape-aware pre-score with small B = 2^PRESELECT_B_BITS
        PreCand pc = buildPreCandWithShape(m, cFloor, c2, residualBits, sieveM, PRESELECT_B_BITS);

        if (heap.size() < PRESELECT_CAP) heap.add(pc);
        else if (pc.preScore < heap.peek().preScore) { heap.poll(); heap.add(pc); }
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
                                          int residualBits, long sieveM, int Bbits) {
        // choose better of floor/ceil c
        BigInteger c, rBest;
        BigInteger c2Floor = cFloor.multiply(cFloor);
        if (c2.equals(c2Floor)) { c = cFloor; rBest = BigInteger.ZERO; }
        else {
            BigInteger rFloor = c2.subtract(c2Floor);
            BigInteger cCeil  = cFloor.add(BigInteger.ONE);
            BigInteger rCeil  = cCeil.multiply(cCeil).subtract(c2);
            if (rCeil.compareTo(rFloor) < 0) { c = cCeil; rBest = rCeil; }
            else { c = cFloor; rBest = rFloor; }
        }

        // VERY shallow factoring for a shape estimate
        List<BigInteger> factors = factorByPrimeDivisionBounded(m, PRE_TRIAL_DIV_LIMIT);
        Split fst = greedySplitToFST(factors, sieveM);

        double pre = eliteEdgeScore(fst.f, fst.s, fst.t, c, sieveM, portionBits, Bbits, residualBits);
        return new PreCand(m, cFloor, c2, residualBits, pre);
    }

    private static final class PreCand {
        final BigInteger m, cFloor, c2;
        final int residualBits;
        final double preScore; // predicted final score (edge-aware)

        private PreCand(BigInteger m, BigInteger cFloor, BigInteger c2, int residualBits, double preScore) {
            this.m = m; this.cFloor = cFloor; this.c2 = c2; this.residualBits = residualBits; this.preScore = preScore;
        }
    }

    // ---------- Stage 2 finalize ----------
    private Poly finalizeOne(PreCand pc, long sieveM) {
        BigInteger cFloor = pc.cFloor;
        BigInteger c2 = pc.c2;

        BigInteger c, rBest;
        BigInteger c2Floor = cFloor.multiply(cFloor);
        if (c2.equals(c2Floor)) { c = cFloor; rBest = BigInteger.ZERO; }
        else {
            BigInteger rFloor = c2.subtract(c2Floor);
            BigInteger cCeil  = cFloor.add(BigInteger.ONE);
            BigInteger rCeil  = cCeil.multiply(cCeil).subtract(c2);
            if (rCeil.compareTo(rFloor) < 0) { c = cCeil; rBest = rCeil; }
            else { c = cFloor; rBest = rFloor; }
        }
        int residualBits = rBest.signum() == 0 ? 0 : rBest.bitLength();

        List<BigInteger> mFactors = factorByPrimeDivision(pc.m);
        Split fst = greedySplitToFST(mFactors, sieveM);
        Poly poly = new Poly(fst.f, fst.s, fst.t, c, residualBits);
        poly.score(sieveM, portionBits);
        return poly;
    }

    private Poly finalizeAndLog(PreCand pc, long sieveM) {
        Poly best = finalizeOne(pc, sieveM);
        log("FSTSearch: best score", String.format("%.2f", best.score), "| fBits", best.f.bitLength(),
                "sBits", best.s.bitLength(), "tBits", best.t.bitLength(),
                "| cBits", best.c.bitLength(), "| linkResidualBits", best.linkResidualBits);
        logEliteSummary(best, sieveM);
        logPredictedEdgeBits(best, sieveM, 20);
        logPredictedEdgeBits(best, sieveM, 30);
        return best;
    }

    // ---------- Elite edge-aware scoring ----------
    /**
     * Predicts edge behavior using slopes and half-window B = 2^Bbits.
     * Uses edge bits (not center) in the RMS + residual penalty.
     */
    private static double eliteEdgeScore(BigInteger f, BigInteger s, BigInteger t, BigInteger c,
                                         long M, int targetBits, int Bbits, int linkResidualBits) {
        if (Bbits < 0) Bbits = 0;
        BigInteger MB = BigInteger.valueOf(M);
        BigInteger fs  = f.multiply(s);
        BigInteger fs2 = fs.multiply(s);

        // Center magnitudes (bits)
        int l1c = s.multiply(MB).subtract(t).abs().bitLength();
        int l2c = s.multiply(MB).add(t).abs().bitLength();
        int l3c = fs2.multiply(MB).subtract(c).abs().bitLength();
        int l4c = fs2.multiply(MB).add(c).abs().bitLength();
        int cons = fs.bitLength() * 2;

        // Slope-driven edge contributions (bits)
        int sBits  = s.bitLength();
        int fsBits = fs.bitLength();

        int l1e = Math.max(l1c, sBits + Bbits);
        int l2e = Math.max(l2c, sBits + Bbits);
        int l3e = Math.max(l3c, 2 * fsBits + Bbits);
        int l4e = Math.max(l4c, 2 * fsBits + Bbits);
        int ce  = cons; // const doesn't change with x

        // RMS vs target at the edge + residual penalty
        double d1 = l1e - targetBits;
        double d2 = l2e - targetBits;
        double d3 = l3e - targetBits;
        double d4 = l4e - targetBits;
        double dc = ce  - targetBits;

        double rms = Math.sqrt((d1*d1 + d2*d2 + d3*d3 + d4*d4 + dc*dc) / 5.0);
        return rms + linkResidualBits * RESIDUAL_WEIGHT;
    }

    // ---------- Behavior logging ----------
    private void logEliteSummary(Poly p, long M) {
        // Slopes and bit sizes
        BigInteger s   = p.s;
        BigInteger fs  = p.f.multiply(p.s);
        BigInteger fs2 = fs.multiply(p.s);

        int sBits  = s.bitLength();
        int fsBits = fs.bitLength();
        int constBits = fsBits * 2;

        log("Elite: sBits=", sBits, " fsBits=", fsBits,
                " | slope(L1/2)=|s|≈2^", sBits,
                " | slope(L3/4)=|f*s^2|≈2^", fs2.bitLength(),
                " | constBits=", constBits);

        // Center and a few specific x values
        logBehaviorAtX(p, M, 0);
        logBehaviorAtX(p, M, 1);
        logBehaviorAtX(p, M, 1L << 10);
        logBehaviorAtX(p, M, 1L << 20);
        logBehaviorAtX(p, M, 1L << 30);
    }

    private void logBehaviorAtX(Poly p, long M, long x) {
        BigInteger MB = BigInteger.valueOf(M);
        BigInteger s  = p.s;
        BigInteger fs2 = p.f.multiply(p.s).multiply(p.s);

        String label = (x == 0) ? "x=0(center)" :
                (x == 1 ? "x=1" :
                        (x == (1L<<10) ? "x=2^10" :
                                (x == (1L<<20) ? "x=2^20" :
                                        (x == (1L<<30) ? "x=2^30" : "x=" + x))));

        // +x
        BigInteger Xp = MB.add(BigInteger.valueOf(x));
        int L1p = s.multiply(Xp).subtract(p.t).abs().bitLength();
        int L2p = s.multiply(Xp).add(p.t).abs().bitLength();
        int L3p = fs2.multiply(Xp).subtract(p.c).abs().bitLength();
        int L4p = fs2.multiply(Xp).add(p.c).abs().bitLength();
        int Cb  = p.f.multiply(p.s).bitLength() * 2;

        if (x == 0) {
            log("Center bits:    L1=", L1p, " L2=", L2p, " L3=", L3p, " L4=", L4p,
                    " | Const=", Cb, " | target=", portionBits);
            return;
        }

        // -x
        BigInteger Xm = MB.subtract(BigInteger.valueOf(x));
        int L1m = s.multiply(Xm).subtract(p.t).abs().bitLength();
        int L2m = s.multiply(Xm).add(p.t).abs().bitLength();
        int L3m = fs2.multiply(Xm).subtract(p.c).abs().bitLength();
        int L4m = fs2.multiply(Xm).add(p.c).abs().bitLength();

        log(label, " (+x) L1=", L1p, " L2=", L2p, " L3=", L3p, " L4=", L4p,
                " | (-x) L1=", L1m, " L2=", L2m, " L3=", L3m, " L4=", L4m,
                " | Const=", Cb, " | target=", portionBits);
    }

    private void logPredictedEdgeBits(Poly p, long M, int Bbits) {
        BigInteger MB = BigInteger.valueOf(M);
        BigInteger s   = p.s;
        BigInteger fs  = p.f.multiply(p.s);
        BigInteger fs2 = fs.multiply(p.s);

        int l1c = s.multiply(MB).subtract(p.t).abs().bitLength();
        int l2c = s.multiply(MB).add(p.t).abs().bitLength();
        int l3c = fs2.multiply(MB).subtract(p.c).abs().bitLength();
        int l4c = fs2.multiply(MB).add(p.c).abs().bitLength();
        int cons = fs.bitLength() * 2;

        int sBits  = s.bitLength();
        int fsBits = fs.bitLength();

        int l1e = Math.max(l1c, sBits + Bbits);
        int l2e = Math.max(l2c, sBits + Bbits);
        int l3e = Math.max(l3c, 2*fsBits + Bbits);
        int l4e = Math.max(l4c, 2*fsBits + Bbits);

        log("Predicted edge bits (M=", M, ", B=2^", Bbits, "):");
        log("  Target:     ", portionBits);
        log("  L1(sx-t):   ", l1e);
        log("  L2(sx+t):   ", l2e);
        log("  L3(fs²x-c): ", l3e);
        log("  L4(fs²x+c): ", l4e);
        log("  Const(f²s²):", cons);
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

            while (true) {
                BigInteger[] qr = n.divideAndRemainder(p);
                if (qr[1].signum() == 0) { res.add(p); n = qr[0]; if (n.equals(BigInteger.ONE)) break; }
                else break;
            }
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
                if (qr[1].signum() == 0) { res.add(p); n = qr[0]; if (n.equals(BigInteger.ONE)) break; }
                else break;
            }
            long p2 = (long) pInt * (long) pInt;
            if (BigInteger.valueOf(p2).compareTo(n) > 0) break;
        }
        if (n.compareTo(BigInteger.ONE) > 0) res.add(n);
        res.sort(Comparator.comparingInt(BigInteger::bitLength).reversed());
        return res;
    }

    // ---------- Greedy F/S/T split ----------
    private Split greedySplitToFST(List<BigInteger> factors, long sieveM) {
        BigInteger f = BigInteger.ONE, s = BigInteger.ONE, t = BigInteger.ONE;
        if (factors.isEmpty()) return new Split(f, s, t);

        for (BigInteger g : factors) {
            double best = Double.POSITIVE_INFINITY;
            int bin = -1;

            double a = scorePreview(f.multiply(g), s, t, sieveM); if (a < best) { best = a; bin = 0; }
            double b = scorePreview(f, s.multiply(g), t, sieveM); if (b < best) { best = b; bin = 1; }
            double c = scorePreview(f, s, t.multiply(g), sieveM); if (c < best) { bin = 2; }

            if (bin == 0) f = f.multiply(g);
            else if (bin == 1) s = s.multiply(g);
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

    // ---------- QS residue wheel ----------
    private static class ScanState {
        final BigInteger m;
        final int residue; // m mod W
        ScanState(BigInteger m, int residue) { this.m = m; this.residue = residue; }
    }

    private static class Wheel {
        final int W;                 // modulus
        final int[] allowedIdx;      // allowed residues in [0..W)
        final int[] nextJump;        // residue r -> distance to next allowed residue
        final int[] prevJump;        // residue r -> distance to prev allowed residue

        Wheel(int W, int[] allowedIdx, int[] nextJump, int[] prevJump) {
            this.W = W; this.allowedIdx = allowedIdx; this.nextJump = nextJump; this.prevJump = prevJump;
        }

        static Wheel build(BigInteger N, int[] primes) {
            final int k = primes.length;
            final int[] fourNmodP = new int[k];
            final boolean[][] isSquare = new boolean[k][];

            for (int i = 0; i < k; i++) {
                int p = primes[i];
                int nmp = N.mod(BigInteger.valueOf(p)).intValue();
                int four = (int) ((4L * nmp) % p);
                fourNmodP[i] = four;

                boolean[] sq = new boolean[p];
                for (int a = 0; a < p; a++) sq[(a * a) % p] = true; // includes 0
                isSquare[i] = sq;
            }

            // Build allowed residues modulo W by incremental CRT
            int W = 1;
            int[] residues = {0};
            for (int i = 0; i < k; i++) {
                int p = primes[i];
                boolean[] sq = isSquare[i];
                int four = fourNmodP[i];

                ArrayList<Integer> okP = new ArrayList<>(p / 2 + 1);
                for (int r = 0; r < p; r++) {
                    int v = r * r - four;
                    v %= p; if (v < 0) v += p;
                    if (sq[v]) okP.add(r);
                }

                int newW = W * p;
                ArrayList<Integer> merged = new ArrayList<>(Math.max(1, residues.length * okP.size() / 2));
                int invWmodP = invMod(W % p, p); // W^{-1} mod p
                for (int a : residues) {
                    int aModP = a % p;
                    for (int r : okP) {
                        int t = (int) (((long) (r - aModP + p) * invWmodP) % p);
                        int x = a + W * t; // unique mod newW
                        merged.add(x);
                    }
                }
                W = newW;
                residues = merged.stream().mapToInt(Integer::intValue).toArray();
            }
            Arrays.sort(residues);

            // Jump tables (circular over [0..W))
            int[] nextJump = new int[W];
            int[] prevJump = new int[W];
            for (int r = 0; r < W; r++) {
                int idx = Arrays.binarySearch(residues, r + 1);
                if (idx < 0) idx = -idx - 1;
                int next = (idx < residues.length) ? residues[idx] : residues[0] + W;
                nextJump[r] = next - r;

                int idxPrev = Arrays.binarySearch(residues, r);
                if (idxPrev >= 0) idxPrev--; else idxPrev = -idxPrev - 2;
                if (idxPrev < 0) idxPrev = residues.length - 1;
                int prev = residues[idxPrev];
                if (prev > r) prev -= W;
                prevJump[r] = r - prev;
            }
            return new Wheel(W, residues, nextJump, prevJump);
        }

        ScanState seedUpFrom(BigInteger m0) {
            int r0 = m0.mod(BigInteger.valueOf(W)).intValue();
            int d = nextJump[r0];
            BigInteger m = (d == 0) ? m0 : m0.add(BigInteger.valueOf(d));
            int r = (r0 + d) % W;
            return new ScanState(m, r);
        }

        ScanState seedDownFrom(BigInteger m0) {
            int r0 = m0.mod(BigInteger.valueOf(W)).intValue();
            int d = prevJump[r0];
            if (m0.signum() <= 0) return null;
            BigInteger m = (d == 0) ? m0 : m0.subtract(BigInteger.valueOf(d));
            if (m.signum() <= 0) return null;
            int r = (r0 - d) % W; if (r < 0) r += W;
            return new ScanState(m, r);
        }

        ScanState stepUp(ScanState s) {
            int d = nextJump[s.residue];
            BigInteger m2 = s.m.add(BigInteger.valueOf(d));
            int r2 = (s.residue + d) % W;
            return new ScanState(m2, r2);
        }

        ScanState stepDown(ScanState s) {
            int d = prevJump[s.residue];
            BigInteger m2 = s.m.subtract(BigInteger.valueOf(d));
            if (m2.signum() <= 0) return null;
            int r2 = (s.residue - d) % W; if (r2 < 0) r2 += W;
            return new ScanState(m2, r2);
        }

        private static int invMod(int a, int p) {
            // a in [0..p), p prime, gcd(a,p)=1
            int t = 0, newT = 1;
            int r = p, newR = a;
            while (newR != 0) {
                int q = r / newR;
                int tmpT = t - q * newT; t = newT; newT = tmpT;
                int tmpR = r - q * newR; r = newR; newR = tmpR;
            }
            if (r != 1) throw new ArithmeticException("invMod: non-invertible");
            if (t < 0) t += p;
            return t;
        }
    }

    // ---------- Utils ----------
    private static int[] sievePrimesUpTo(int limit) {
        boolean[] isComposite = new boolean[limit + 1];
        for (int i = 2; i * i <= limit; i++) {
            if (!isComposite[i]) for (int j = i * i; j <= limit; j += i) isComposite[j] = true;
        }
        int count = 0;
        for (int i = 2; i <= limit; i++) if (!isComposite[i]) count++;
        int[] primes = new int[count];
        int idx = 0;
        for (int i = 2; i <= limit; i++) if (!isComposite[i]) primes[idx++] = i;
        return primes;
    }

    // ---------- Result type ----------
    public static class Poly {
        public final BigInteger f, s, t, c;
        public final int linkResidualBits; // 0 for exact; >0 for approximate
        public double score;

        public Poly(BigInteger f, BigInteger s, BigInteger t, BigInteger c, int linkResidualBits) {
            this.f = f; this.s = s; this.t = t; this.c = c; this.linkResidualBits = linkResidualBits;
        }

        private static double sq(double x) { return x * x; }

        /**
         * Score k=5 portions:
         * L1 = |s*M - t|,  L2 = |s*M + t|,
         * L3 = |f*s*s*M - c|,  L4 = |f*s*s*M + c|,
         * const = (f*s)^2.
         * Minimizes RMS vs 'targetBits', plus residual penalty.
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

        @Override public String toString() {
            return "Poly{f=" + f + ", s=" + s + ", t=" + t + ", c=" + c +
                    ", linkResidualBits=" + linkResidualBits + ", score=" + score + "}";
        }
    }
}
