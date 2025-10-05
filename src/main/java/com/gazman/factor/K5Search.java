package com.gazman.factor;

import com.gazman.math.MathUtils;
import com.gazman.math.SqrRoot;

import java.math.BigInteger;
import java.util.*;

/**
 * Fixed-time search for (f, s, t, c) tuples for k=5 with PROGRAMMATIC bit-budgeting:
 *  - Determine s,t,f bit windows from: xBits, portion target P=|N|/5, and the constraint
 *        cBits ≈ 2*sumBits  and  c ≈ f*s^2*x_center .
 *  - Build s,t near those bit targets; f approximates value ~ floor(sqrt(2N))/(s*t).
 *  - Hensel/CRT for c ≡ sqrt(2N) (mod q), with q=(fst)^2 (centered representative).
 *  - Accept only tuples for which the implied x-centers x_t=t/s and x_c=c/(f s^2) both land
 *    in the QS x-window and predict L1,L3 sizes near the portion band.
 *  - Telemetry logs targets, windows, and best tuple summary.
 *  - NEW: per-second K5SearchSpeed logs (attempts/s, hensel/s, accepted/s, keep size, best period stats).
 */
public class K5Search extends Logger {

    // ---- Public result tuple ----
    public static final class Tuple {
        public final BigInteger f, s, t, c, q, db;
        public final int fBits, sBits, tBits, cBits, qBits, dbBits;
        public final long xT, xC;        // centers (rounded) for L1/L3
        public final int l1BitsPred, l3BitsPred;
        public final int primesCount;    // distinct primes in f∪s∪t
        public final int sumBits;        // fBits + sBits + tBits

        public Tuple(BigInteger f, BigInteger s, BigInteger t, BigInteger c,
                     BigInteger q, BigInteger db,
                     int fBits, int sBits, int tBits, int cBits, int qBits, int dbBits,
                     long xT, long xC, int l1BitsPred, int l3BitsPred,
                     int primesCount, int sumBits) {
            this.f = f; this.s = s; this.t = t; this.c = c;
            this.q = q; this.db = db;
            this.fBits = fBits; this.sBits = sBits; this.tBits = tBits;
            this.cBits = cBits; this.qBits = qBits; this.dbBits = dbBits;
            this.xT = xT; this.xC = xC;
            this.l1BitsPred = l1BitsPred; this.l3BitsPred = l3BitsPred;
            this.primesCount = primesCount; this.sumBits = sumBits;
        }
    }

    // ---- Internal helper: product with exact prime list ----
    private static final class Product {
        final BigInteger val;
        final ArrayList<Integer> primes; // distinct primes used (from pool)
        final int bits;
        Product(BigInteger val, ArrayList<Integer> primes) {
            this.val = val;
            this.primes = primes;
            this.bits = val.bitLength();
        }
    }

    // ---- Telemetry ----
    private final class Stats {
        int poolSize;
        int nbits;
        int highestPrime;
        int xBits;

        int portionBits;          // ≈ |N|/5
        int fs2BitsTarget;        // target bits for f*s^2
        int cBitsTarget;          // fs2BitsTarget + xBits
        int sumBitsTarget;        // ≈ (cBitsTarget + 1)/2

        int sTarget, tTarget, fTargetBits;
        int sMin, sMax, tMin, tMax;
        int minFBits;             // lower bound for f’s bits (stability)
        int fTolAbs;              // tolerance for f bits around target

        long loops;
        long sBuilt, sFailBits;
        long fBuilt, fFailBits;
        long tBuilt, tFailBits;
        long tuplesTried;
        long henselOk, henselFail;
        long gapNotDivisible;
        long accepted;

        int bestDbBits = Integer.MAX_VALUE;
        int bestL1Pred = Integer.MAX_VALUE;
        int bestL3Pred = Integer.MAX_VALUE;

        // ---- per-second speed accounting ----
        long lastSpeedLogMs;
        long attemptsPeriod, henselPeriod, acceptedPeriod;
        int bestDbPeriod = Integer.MAX_VALUE;
        int bestL1Period = Integer.MAX_VALUE;
        int bestL3Period = Integer.MAX_VALUE;

        void header(List<Integer> preview) {
            log("K5Search: pool size", poolSize, "| first primes", preview);
            log("K5Search: Nbits", nbits,
                    "xBits", xBits, "portion", portionBits,
                    "fs2Bits*", fs2BitsTarget, "cBits*", cBitsTarget, "sumBits*", sumBitsTarget);
            log("K5Search: targets | sBits~", sTarget, "tBits~", tTarget, "fBits~", fTargetBits,
                    " | s in [", sMin, ",", sMax, "] t in [", tMin, ",", tMax, "]",
                    " minFBits", minFBits, "fTolAbs", fTolAbs);
        }

        void maybeSpeedLog(PriorityQueue<Tuple> keep) {
            long now = System.currentTimeMillis();
            if (lastSpeedLogMs == 0) {
                lastSpeedLogMs = now;
                return;
            }
            long dt = now - lastSpeedLogMs;
            if (dt < 1000) return;

            double sec = dt / 1000.0;
            long attPerS = Math.round(attemptsPeriod / sec);
            long henPerS = Math.round(henselPeriod / sec);
            long accPerS = Math.round(acceptedPeriod / sec);

            int dbBest = bestDbPeriod == Integer.MAX_VALUE ? -1 : bestDbPeriod;
            int l1Best = bestL1Period == Integer.MAX_VALUE ? -1 : bestL1Period;
            int l3Best = bestL3Period == Integer.MAX_VALUE ? -1 : bestL3Period;

            log("K5SearchSpeed:",
                    attPerS / 1000, "kAtt/s,",
                    "hensel", henPerS / 1000, "k/s,",
                    "accept", accPerS, "/s |",
                    "keep", keep.size(),
                    "| best(dbBits,L1~,L3~)", dbBest, l1Best, l3Best);

            // reset period counters
            attemptsPeriod = henselPeriod = acceptedPeriod = 0;
            bestDbPeriod = bestL1Period = bestL3Period = Integer.MAX_VALUE;
            lastSpeedLogMs = now;
        }

        void footer(long dt, PriorityQueue<Tuple> keep) {
            if (keep.isEmpty()) {
                log("K5Search: no tuples found — telemetry");
                log("  loops", loops,
                        "| sBuilt", sBuilt, "sFailBits", sFailBits,
                        "| fBuilt", fBuilt, "fFailBits", fFailBits,
                        "| tBuilt", tBuilt, "tFailBits", tFailBits);
                log("  tuplesTried", tuplesTried,
                        "| henselOk", henselOk, "henselFail", henselFail,
                        "| gapNotDivisible", gapNotDivisible,
                        "| bestDbBits", (bestDbBits==Integer.MAX_VALUE? -1:bestDbBits),
                        "| bestL1Pred", (bestL1Pred==Integer.MAX_VALUE? -1:bestL1Pred),
                        "| bestL3Pred", (bestL3Pred==Integer.MAX_VALUE? -1:bestL3Pred));
            } else {
                Tuple b = keep.stream()
                        .min(Comparator
                                .comparingInt((Tuple t)-> Math.abs(t.l1BitsPred - portionBits)
                                        + Math.abs(t.l3BitsPred - portionBits))
                                .thenComparingInt(t -> t.dbBits))
                        .get();
                log("K5Search: time", dt, "ms | attempts", tuplesTried,
                        "| henselOk", henselOk, "| accepted", accepted,
                        "| best tuple -> s", b.sBits, "t", b.tBits, "f", b.fBits,
                        "| c", b.cBits, "q", b.qBits, "db", b.dbBits,
                        "| xT", b.xT, "xC", b.xC,
                        "| L1~", b.l1BitsPred, "L3~", b.l3BitsPred);
            }
        }
    }

    private final BigInteger N;
    private final BigInteger N2;
    private final int[] primeBaseInt;
    private final Random rnd = new Random(1337);

    // acceptance knobs
    private static final int BAND_LO = 8;     // portion - BAND_LO
    private static final int BAND_HI = 10;    // portion + BAND_HI
    private static final int DB_BITS_MAX = 140; // loose cap on dbBits

    public K5Search(BigInteger N, int[] primeBaseInt) {
        this.N = N;
        this.N2 = N.shiftLeft(1);
        this.primeBaseInt = primeBaseInt;
    }

    /** Run for a fixed ms, return up to topK tuples. Also logs telemetry if empty. */
    public List<Tuple> runFixedTime(int msBudget, int topK) {
        long start = System.currentTimeMillis();
        long end = start + msBudget;

        // pool of odd primes p with (2N|p)=+1 (ceiling ~2e6)
        ArrayList<Integer> good = new ArrayList<>(8192);
        int highest = 3;
        for (int p : primeBaseInt) {
            if (p <= 2) continue;
            highest = Math.max(highest, p);
            if (p > 2_000_000) break;
            if (MathUtils.isRootInQuadraticResidues(N2, BigInteger.valueOf(p))) {
                good.add(p);
            }
        }

        Stats st = new Stats();
        st.poolSize = good.size();
        st.nbits = N.bitLength();
        st.highestPrime = highest;
        st.xBits = Integer.toBinaryString(highest).length() - 1; // floor(log2(highest))

        if (good.isEmpty()) {
            log("K5Search: pool empty — no primes with (2N|p)=+1 under 2e6.");
            return List.of();
        } else {
            // ---- PROGRAMMATIC BIT BUDGET ----
            st.portionBits   = Math.max(24, st.nbits / 5);
            st.fs2BitsTarget = Math.max(12, st.portionBits - st.xBits);
            st.cBitsTarget   = st.fs2BitsTarget + st.xBits;
            st.sumBitsTarget = (st.cBitsTarget + 1) / 2;

            // choose s,t,f targets: enforce tBits ≈ sBits + xBits; f fills the rest
            st.minFBits = Math.max(10, st.nbits / 16);
            st.fTolAbs  = 6;

            // pick s so that fBits stays above min; sTarget from splitting leftover
            int sTarget = Math.max(8, (st.sumBitsTarget - st.xBits - st.minFBits) / 2);
            int fBits   = st.sumBitsTarget - (sTarget + (sTarget + st.xBits));
            if (fBits < st.minFBits) {
                // push s down to make room for f
                int deficit = st.minFBits - fBits;
                sTarget = Math.max(8, sTarget - (deficit + 1) / 2);
                fBits = st.sumBitsTarget - (sTarget + (sTarget + st.xBits));
            }
            st.sTarget = sTarget;
            st.tTarget = sTarget + st.xBits;
            st.fTargetBits = Math.max(st.minFBits, fBits);

            int slackS = Math.max(2, st.sTarget / 5);
            int slackT = Math.max(2, st.tTarget / 5);
            st.sMin = Math.max(8, st.sTarget - slackS);
            st.sMax = st.sTarget + slackS;
            st.tMin = Math.max(8, st.tTarget - slackT);
            st.tMax = st.tTarget + slackT;

            List<Integer> preview = good.subList(0, Math.min(10, good.size()));
            st.header(preview);

            PriorityQueue<Tuple> keep = new PriorityQueue<>(Comparator
                    .comparingInt((Tuple t) ->
                            Math.abs(t.l1BitsPred - st.portionBits)
                                    + Math.abs(t.l3BitsPred - st.portionBits))
                    .thenComparingInt(t -> t.dbBits)
                    .reversed()); // worst on top

            Tuple bestAny = null;
            HashSet<Integer> usedSet = new HashSet<>(64);
            BigInteger cApprox = SqrRoot.bigIntSqRootFloor(N2);

            // init speed timer
            st.lastSpeedLogMs = System.currentTimeMillis();

            while (System.currentTimeMillis() < end) {
                st.loops++;
                usedSet.clear();

                // s near target
                Product sP = growFactorBits(good, usedSet, st.sMin, st.sMax);
                if (sP.val.equals(BigInteger.ONE)) { st.sFailBits++; continue; }
                st.sBuilt++;

                // t near s + xBits
                Product tP = growFactorBits(good, usedSet, st.tMin, st.tMax);
                if (tP.val.equals(BigInteger.ONE)) { st.tFailBits++; continue; }
                st.tBuilt++;

                // f by VALUE: target ~ floor(sqrt(2N)) / (s * t)
                BigInteger fValTarget = cApprox.divide(sP.val.multiply(tP.val));
                if (fValTarget.signum() <= 0) fValTarget = BigInteger.ONE.shiftLeft(st.fTargetBits);
                Product fP = growFactorApproxValue(good, usedSet, fValTarget, st.minFBits, st.fTolAbs);
                if (fP.val.equals(BigInteger.ONE)) { st.fFailBits++; continue; }
                st.fBuilt++;

                st.tuplesTried++;
                st.attemptsPeriod++;

                BigInteger sBI = sP.val, tBI = tP.val, fBI = fP.val;
                int sumBits = fP.bits + sP.bits + tP.bits;

                // q = (fst)^2
                BigInteger q = fBI.multiply(sBI).multiply(tBI);
                q = q.multiply(q);

                // distinct primes used
                ArrayList<Integer> primes = new ArrayList<>(sP.primes.size() + fP.primes.size() + tP.primes.size());
                primes.addAll(sP.primes); primes.addAll(fP.primes); primes.addAll(tP.primes);
                primes = new ArrayList<>(new LinkedHashSet<>(primes));

                long xT = safeRoundDiv(tBI, sBI);                   // ≈ t/s
                BigInteger fs2 = fBI.multiply(sBI).multiply(sBI);   // f*s^2
                BigInteger targetC = fs2.multiply(BigInteger.valueOf(xT));

                // plain Hensel/CRT (center to 0)
                BigInteger c0 = sqrtN2ModCompositeSquareRecenter(q, primes, targetC);
                if (c0 == null) { st.henselFail++; st.maybeSpeedLog(keep); continue; }
                st.henselOk++; st.henselPeriod++;

                // db = |2N - c^2| / q
                BigInteger c2 = c0.multiply(c0);
                BigInteger gap = N2.subtract(c2).abs();
                if (!gap.mod(q).equals(BigInteger.ZERO)) { st.gapNotDivisible++; st.maybeSpeedLog(keep); continue; }
                BigInteger db = peelTiny(gap.divide(q), 512);

                int dbBits   = db.bitLength();
                int cBits    = c0.bitLength();

                // x centers (rounded)
                long xC = safeRoundDiv(c0.abs(), fs2);

                // discard if centers out of QS window
                if (xT < 0 || xT >= (long) st.highestPrime) { st.maybeSpeedLog(keep); continue; }
                if (xC < 0 || xC >= (long) st.highestPrime) { st.maybeSpeedLog(keep); continue; }

                // predicted sizes at x ≈ xT
                long dx = Math.abs(xC - xT);
                if (dx == 0) dx = 1; // avoid zero; we still need some width
                BigInteger l1Pred = sBI.multiply(BigInteger.valueOf(dx));       // |s*x - t| ~ s*|Δx|
                BigInteger l3Pred = fs2.multiply(BigInteger.valueOf(dx));       // |fs^2*x - c| ~ fs^2*|Δx|
                int l1BitsPred = l1Pred.bitLength();
                int l3BitsPred = l3Pred.bitLength();

                st.bestDbBits = Math.min(st.bestDbBits, dbBits);
                st.bestL1Pred = Math.min(st.bestL1Pred, l1BitsPred);
                st.bestL3Pred = Math.min(st.bestL3Pred, l3BitsPred);

                // period bests
                st.bestDbPeriod = Math.min(st.bestDbPeriod, dbBits);
                st.bestL1Period = Math.min(st.bestL1Period, l1BitsPred);
                st.bestL3Period = Math.min(st.bestL3Period, l3BitsPred);

                // portion band check
                boolean l1Ok = (l1BitsPred >= st.portionBits - BAND_LO) && (l1BitsPred <= st.portionBits + BAND_HI);
                boolean l3Ok = (l3BitsPred >= st.portionBits - BAND_LO) && (l3BitsPred <= st.portionBits + BAND_HI);

                // db sanity (don’t waste on huge leftovers)
                if (dbBits > DB_BITS_MAX && rnd.nextInt(8) != 0) { st.maybeSpeedLog(keep); continue; }

                Tuple tup = new Tuple(
                        fP.val, sP.val, tP.val, c0, q, db,
                        fP.bits, sP.bits, tP.bits, cBits, q.bitLength(), dbBits,
                        xT, xC, l1BitsPred, l3BitsPred,
                        primes.size(), sumBits
                );

                if (bestAny == null || worseThan(bestAny, tup, st.portionBits)) {
                    bestAny = tup;
                }

                if (l1Ok && l3Ok) {
                    st.accepted++; st.acceptedPeriod++;
                    if (keep.size() < topK) {
                        keep.add(tup);
                    } else {
                        Tuple worst = keep.peek();
                        if (!worseThan(tup, worst, st.portionBits)) { keep.poll(); keep.add(tup); }
                    }
                }

                st.maybeSpeedLog(keep);
            }

            long dt = System.currentTimeMillis() - start;

            ArrayList<Tuple> out = new ArrayList<>(keep);
            out.sort(Comparator
                    .comparingInt((Tuple t) -> Math.abs(t.l1BitsPred - st.portionBits)
                            + Math.abs(t.l3BitsPred - st.portionBits))
                    .thenComparingInt(t -> t.dbBits));

            st.footer(dt, keep);

            if (out.isEmpty()) {
                if (bestAny != null) {
                    log("SearchSummary[K=5]: f", bestAny.f, "(" + bestAny.fBits + ")",
                            "s", bestAny.s, "(" + bestAny.sBits + ")",
                            "t", bestAny.t, "(" + bestAny.tBits + ")",
                            "c", bestAny.c, "(" + bestAny.cBits + ")",
                            "db", bestAny.db, "(" + bestAny.dbBits + ")",
                            "xT", bestAny.xT, "xC", bestAny.xC,
                            "L1~", bestAny.l1BitsPred, "L3~", bestAny.l3BitsPred);
                }
            } else {
                log("K5Search: top", Math.min(topK, out.size()), "tuples (close to portion band)");
                for (int i = 0; i < Math.min(topK, out.size()); i++) {
                    Tuple t = out.get(i);
                    log("  #" + (i+1) + ": f", t.fBits, "bits, s", t.sBits, "bits, t", t.tBits, "bits,",
                            "c", t.cBits, "bits | Q", t.qBits, "bits | db", t.dbBits, "bits |",
                            "xT", t.xT, "xC", t.xC, "| L1~", t.l1BitsPred, "L3~", t.l3BitsPred,
                            " | primes=", t.primesCount);
                }
            }
            return out;
        }
    }

    private static boolean worseThan(Tuple a, Tuple b, int portion) {
        int as = Math.abs(a.l1BitsPred - portion) + Math.abs(a.l3BitsPred - portion);
        int bs = Math.abs(b.l1BitsPred - portion) + Math.abs(b.l3BitsPred - portion);
        if (as != bs) return as > bs;
        return a.dbBits > b.dbBits;
    }

    /** Grow a factor to hit a bit window. Distinct primes (no reuse with usedSet). */
    private Product growFactorBits(ArrayList<Integer> pool, HashSet<Integer> usedSet,
                                   int minBits, int maxBits) {
        if (maxBits < 1) return new Product(BigInteger.ONE, new ArrayList<>());
        BigInteger x = BigInteger.ONE;
        ArrayList<Integer> primes = new ArrayList<>(16);

        for (int tries = 0; tries < 256 && x.bitLength() < minBits; tries++) {
            int idx = new Random().nextInt(pool.size());
            int p = pool.get(idx);
            if (usedSet.contains(p)) continue;
            BigInteger nx = x.multiply(BigInteger.valueOf(p));
            if (nx.bitLength() > maxBits) continue;
            x = nx;
            primes.add(p);
            usedSet.add(p);
        }

        if (x.bitLength() < minBits) {
            return new Product(BigInteger.ONE, new ArrayList<>());
        }
        return new Product(x, primes);
    }

    /**
     * Build a product approximating 'target' by VALUE (not only bits),
     * using distinct primes from 'pool' not in usedSet.
     * Accept only when candidate’s bits lie in [minBits .. target±fTolAbs].
     * If no candidate lands in window, returns ONE (fail attempt).
     */
    private Product growFactorApproxValue(ArrayList<Integer> pool,
                                          HashSet<Integer> usedSet,
                                          BigInteger target, int minBits, int fTolAbs) {
        if (target.signum() <= 0) return new Product(BigInteger.ONE, new ArrayList<>());

        final int targetBits = Math.max(minBits, target.bitLength());
        final int minAcceptBits = Math.max(minBits, targetBits - fTolAbs);
        final int maxAcceptBits = targetBits + fTolAbs;

        record Cand(BigInteger val, int bits, ArrayList<Integer> primes, int scoreBits, int deltaLen) {}

        ArrayList<Cand> beam = new ArrayList<>(8);
        beam.add(new Cand(BigInteger.ONE, 0, new ArrayList<>(), targetBits, target.bitLength())); // seed

        final int MAX_ROUNDS = 10;
        final int LOCAL_WIDTH = 10;
        int poolLen = pool.size();
        Random rnd = new Random(1337);

        Cand bestInWindow = null;

        for (int round = 0; round < MAX_ROUNDS; round++) {
            ArrayList<Cand> next = new ArrayList<>(LOCAL_WIDTH * 4);

            for (Cand c : beam) {
                for (int t = 0; t < 16; t++) {
                    int idx = rnd.nextInt(poolLen);
                    int p = pool.get(idx);
                    if (usedSet.contains(p)) continue;
                    if (c.primes.contains(p)) continue;
                    BigInteger newVal = c.val.multiply(BigInteger.valueOf(p));
                    int nb = newVal.bitLength();

                    if (nb > maxAcceptBits + 2) continue;

                    int scoreBits = Math.abs(nb - targetBits);
                    int deltaLen = (nb == 0) ? targetBits : newVal.subtract(target).abs().bitLength();

                    Cand cand = new Cand(newVal, nb, withAdded(c.primes, p), scoreBits, deltaLen);
                    next.add(cand);

                    if (nb >= minAcceptBits && nb <= maxAcceptBits) {
                        if (bestInWindow == null
                                || cand.deltaLen < bestInWindow.deltaLen
                                || (cand.deltaLen == bestInWindow.deltaLen && cand.scoreBits < bestInWindow.scoreBits)) {
                            bestInWindow = cand;
                        }
                    }
                }
            }

            if (bestInWindow != null) {
                for (int p : bestInWindow.primes) usedSet.add(p);
                return new Product(bestInWindow.val, new ArrayList<>(bestInWindow.primes));
            }

            if (next.isEmpty()) break;

            next.sort(Comparator
                    .comparingInt((Cand c) -> c.scoreBits)
                    .thenComparingInt(c -> c.deltaLen));
            beam = new ArrayList<>(next.subList(0, Math.min(LOCAL_WIDTH, next.size())));
        }

        return new Product(BigInteger.ONE, new ArrayList<>());
    }

    private static ArrayList<Integer> withAdded(ArrayList<Integer> src, int p) {
        ArrayList<Integer> dst = new ArrayList<>(src.size() + 1);
        dst.addAll(src);
        dst.add(p);
        return dst;
    }

    /** Light preview peel by the tiniest base primes. */
    private BigInteger peelTiny(BigInteger n, int howMany) {
        BigInteger x = n;
        int taken = 0;
        for (int p : primeBaseInt) {
            if (p < 2) continue;
            if (taken >= howMany) break;
            BigInteger P = BigInteger.valueOf(p);
            while (x.mod(P).equals(BigInteger.ZERO)) x = x.divide(P);
            taken++;
            if (x.equals(BigInteger.ONE)) return x;
        }
        return x;
    }

    // ---------- CRT/Hensel (centered to 0) ----------
    // ---------- CRT/Hensel (recenter to an arbitrary target) ----------
    /**
     * Returns c ≡ sqrt(2N) (mod q) with q = ∏ p_i^2 (composite square),
     * choosing the representative closest to 'target' (instead of 0).
     * 'primes' must contain the distinct odd primes used in (f,s,t).
     */
    private BigInteger sqrtN2ModCompositeSquareRecenter(BigInteger q,
                                                        ArrayList<Integer> primes,
                                                        BigInteger target) {
        BigInteger r = BigInteger.ZERO; // residue mod m
        BigInteger m = BigInteger.ONE;

        for (int p : primes) {
            if (p <= 2) return null;
            BigInteger P = BigInteger.valueOf(p);
            BigInteger p2 = P.multiply(P);
            BigInteger a = N.shiftLeft(1).mod(p2);  // 2N mod p^2

            // sqrt mod p
            long[] roots = MathUtils.ressol(p, a.mod(P).longValue());
            long r0 = -1;
            for (long v : roots) { if (v >= 0) { r0 = v; break; } }
            if (r0 < 0) return null;
            BigInteger s = BigInteger.valueOf(r0);

            // Hensel to p^2:
            BigInteger s2 = s.multiply(s).mod(p2);
            BigInteger diff = a.subtract(s2);
            if (diff.signum() < 0) diff = diff.add(p2);
            BigInteger rhs = diff.divide(P).mod(P);
            BigInteger inv2s = s.shiftLeft(1).mod(P).modInverse(P);
            BigInteger k = rhs.multiply(inv2s).mod(P);
            BigInteger sLift = s.add(k.multiply(P)).mod(p2);    // sqrt mod p^2
            BigInteger sLiftNeg = p2.subtract(sLift).mod(p2);   // the other sign

            // combine each sign with current residue and pick closer to 'target'
            BigInteger rResidue = r.mod(m);
            BigInteger r1 = crtCombine(rResidue, m, sLift, p2);     if (r1 == null) return null;
            BigInteger r2 = crtCombine(rResidue, m, sLiftNeg, p2);  if (r2 == null) return null;

            BigInteger mNew = m.multiply(p2);

            BigInteger c1 = recenterToTarget(r1, mNew, target);
            BigInteger c2 = recenterToTarget(r2, mNew, target);

            // pick the representative closer to 'target'
            BigInteger d1 = c1.subtract(target).abs();
            BigInteger d2 = c2.subtract(target).abs();
            r = (d1.compareTo(d2) <= 0) ? c1 : c2;
            m = mNew;
        }

        if (!m.equals(q)) return null;
        return r;
    }


    /** Combine (r mod m) and (x mod n) with gcd=1; return null if not coprime. */
    private static BigInteger crtCombine(BigInteger r, BigInteger m, BigInteger x, BigInteger n) {
        BigInteger g = m.gcd(n);
        if (!g.equals(BigInteger.ONE)) return null;
        BigInteger inv = m.modInverse(n);
        BigInteger t = x.subtract(r).mod(n);
        BigInteger coeff = t.multiply(inv).mod(n);
        return r.add(m.multiply(coeff));
    }

    /** Re-center residue 'r (mod M)' to the nearest integer to 'target'. */
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

    private static long safeRoundDiv(BigInteger num, BigInteger den) {
        if (den.signum() == 0) return Long.MAX_VALUE;
        BigInteger[] qr = num.add(den.shiftRight(1)).divideAndRemainder(den);
        BigInteger q = qr[0];
        if (q.compareTo(BigInteger.valueOf(Long.MAX_VALUE)) > 0) return Long.MAX_VALUE;
        if (q.compareTo(BigInteger.ZERO) < 0) return Long.MIN_VALUE;
        return q.longValue();
    }
}
