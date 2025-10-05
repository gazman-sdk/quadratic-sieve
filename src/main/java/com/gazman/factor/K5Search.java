package com.gazman.factor;

import com.gazman.math.MathUtils;
import com.gazman.math.SqrRoot;

import java.math.BigInteger;
import java.util.*;

/**
 * Fixed-time search for (f, s, t, c) tuples for k=5 with:
 *  - bit-budget targeting (~|N|/5 per portion),
 *  - recentering c near sqrt(2N) to minimize |2N - c^2|,
 *  - alignment scoring so t/s ≈ c/(f s^2),
 *  - preference for many distinct small primes in f∪s∪t,
 *  - light peel of db by the tiniest primes for preview.
 *
 * Outputs top tuples (by dbBits, then alignment, then distinct-prime count).
 * Does not alter sieving (report-only).
 */
public class K5Search extends Logger {

    public static final class Tuple {
        public final BigInteger f, s, t, c, q, db;
        public final int fBits, sBits, tBits, cBits, qBits, dbBits;
        public final int alignBits;   // log2(|t f s^2 - c s|) - log2(f s^3)  (smaller is better)
        public final int primesCount; // distinct primes in f∪s∪t

        public Tuple(BigInteger f, BigInteger s, BigInteger t, BigInteger c,
                     BigInteger q, BigInteger db,
                     int fBits, int sBits, int tBits, int cBits, int qBits, int dbBits,
                     int alignBits, int primesCount) {
            this.f = f; this.s = s; this.t = t; this.c = c;
            this.q = q; this.db = db;
            this.fBits = fBits; this.sBits = sBits; this.tBits = tBits;
            this.cBits = cBits; this.qBits = qBits; this.dbBits = dbBits;
            this.alignBits = alignBits; this.primesCount = primesCount;
        }
    }

    private final BigInteger N;
    private final BigInteger N2;
    private final int[] primeBaseInt;
    private final Random rnd = new Random(1337);

    public K5Search(BigInteger N, int[] primeBaseInt) {
        this.N = N;
        this.N2 = N.shiftLeft(1);
        this.primeBaseInt = primeBaseInt;
    }

    /** Run for a fixed ms, return up to topK tuples. */
    public List<Tuple> runFixedTime(int msBudget, int topK) {
        long end = System.currentTimeMillis() + msBudget;

        // pool of odd primes p with (2N|p)=+1
        ArrayList<Integer> good = new ArrayList<>(4096);
        for (int p : primeBaseInt) {
            if (p <= 2) continue;
            if (p > 200_000) break; // keep small for cheap lifting
            if (MathUtils.isRootInQuadraticResidues(N2, BigInteger.valueOf(p))) {
                good.add(p);
            }
        }
        if (good.isEmpty()) return List.of();

        final int Nbits = N.bitLength();
        final int portion = Math.max(8, Nbits / 5);   // target |N|/5 bits

        // Windows for s, f, t bits (aim N/5; f satisfies f+2s ≈ N/5)
        final int sTarget = portion;
        final int sTol = Math.max(4, sTarget / 6);
        final int sMin = Math.max(4, sTarget - sTol);
        final int sMax = sTarget + sTol;

        // We'll pick s first, then set fTarget = max(0, portion - 2*sBits) ± small tol
        final int fTolAbs = 4;

        // t: keep similar scale to s
        final int tMin = sMin, tMax = sMax;

        // beam of best tuples
        PriorityQueue<Tuple> best = new PriorityQueue<>(Comparator
                .comparingInt((Tuple t) -> t.dbBits)
                .thenComparingInt(t -> t.alignBits)
                .thenComparingInt(t -> -t.primesCount)
                .thenComparingInt(t -> t.cBits)
                .reversed()); // so we can drop worst cheaply

        // work buffers
        ArrayList<Integer> used = new ArrayList<>(64);
        HashSet<Integer> usedSet = new HashSet<>(64);

        while (System.currentTimeMillis() < end) {
            used.clear();
            usedSet.clear();

            // choose s near target
            BigInteger s = growFactor(good, used, usedSet, sMin, sMax);
            if (s.equals(BigInteger.ONE)) continue;
            int sBits = s.bitLength();

            // choose f so fBits ≈ portion - 2*sBits (clamped)
            int fTarget = Math.max(0, portion - 2 * sBits);
            int fMin = Math.max(0, fTarget - fTolAbs);
            int fMax = fTarget + fTolAbs;
            BigInteger f = growFactor(good, used, usedSet, fMin, fMax);
            if (f.equals(BigInteger.ONE)) continue;

            // choose t near s scale
            BigInteger t = growFactor(good, used, usedSet, tMin, tMax);
            if (t.equals(BigInteger.ONE)) continue;

            int distinctPrimes = usedSet.size();

            BigInteger q = f.multiply(s).multiply(t);
            q = q.multiply(q); // (fst)^2

            // sqrt(2N) mod q via Tonelli/Hensel to p^2 and CRT
            BigInteger c0 = sqrtN2ModCompositeSquare(q, primesOf(f, s, t));
            if (c0 == null) continue;

            // recenter c near sqrt(2N)
            BigInteger c = recenterC(c0, q, N2);

            // compute db = |2N - c^2| / q (exactly divisible)
            BigInteger c2 = c.multiply(c);
            BigInteger gap = N2.subtract(c2).abs();
            if (!gap.mod(q).equals(BigInteger.ZERO)) continue; // safety
            BigInteger db = gap.divide(q);

            // quick peel by the tiniest primes as a preview
            db = peelTiny(db, 512);

            // alignment: Δx ≈ |t/s - c/(f s^2)|; compare bitwise numerator vs denominator
            BigInteger alignNum = t.multiply(f).multiply(s).multiply(s).subtract(c.multiply(s)).abs(); // |t f s^2 - c s|
            BigInteger alignDen = f.multiply(s).multiply(s).multiply(s);                               // f s^3
            int alignBits = alignNum.equals(BigInteger.ZERO)
                    ? Integer.MIN_VALUE
                    : (alignNum.bitLength() - alignDen.bitLength());

            Tuple tup = new Tuple(
                    f, s, t, c, q, db,
                    f.bitLength(), s.bitLength(), t.bitLength(), c.bitLength(), q.bitLength(), db.bitLength(),
                    alignBits, distinctPrimes
            );

            if (best.size() < topK) {
                best.add(tup);
            } else {
                Tuple worst = best.peek();
                if (worseThan(tup, worst)) {
                    // keep
                } else {
                    best.poll();
                    best.add(tup);
                }
            }
        }

        ArrayList<Tuple> out = new ArrayList<>(best);
        out.sort(Comparator
                .comparingInt((Tuple t) -> t.dbBits)
                .thenComparingInt(t -> t.alignBits)
                .thenComparingInt(t -> -t.primesCount)
                .thenComparingInt(t -> t.cBits));
        return out;
    }

    private static boolean worseThan(Tuple a, Tuple b) {
        if (a.dbBits != b.dbBits) return a.dbBits > b.dbBits;
        if (a.alignBits != b.alignBits) return a.alignBits > b.alignBits;
        if (a.primesCount != b.primesCount) return a.primesCount < b.primesCount;
        return a.cBits > b.cBits;
    }

    /** Grow a factor as a product of distinct primes so its bit length hits [minBits, maxBits]. */
    private BigInteger growFactor(ArrayList<Integer> pool, ArrayList<Integer> used, HashSet<Integer> usedSet,
                                  int minBits, int maxBits) {
        if (maxBits < 1) return BigInteger.ONE;
        BigInteger x = BigInteger.ONE;
        int tries = 0;
        while (x.bitLength() < minBits && tries < 128) {
            int idx = rnd.nextInt(pool.size());
            int p = pool.get(idx);
            if (usedSet.contains(p)) { tries++; continue; }
            BigInteger nx = x.multiply(BigInteger.valueOf(p));
            if (nx.bitLength() > maxBits) { tries++; continue; }
            x = nx;
            used.add(p);
            usedSet.add(p);
        }
        return x;
    }

    /** Distinct primes used by f, s, t (we built them from small primes). */
    private static ArrayList<Integer> primesOf(BigInteger f, BigInteger s, BigInteger t) {
        ArrayList<Integer> ps = new ArrayList<>(64);
        factorIntoDistinct(f, ps);
        factorIntoDistinct(s, ps);
        factorIntoDistinct(t, ps);
        return ps;
    }

    private static void factorIntoDistinct(BigInteger n, ArrayList<Integer> out) {
        BigInteger x = n;
        for (int p = 3; p <= 1_000_000 && (long)p * (long)p <= x.longValue(); p += 2) {
            BigInteger P = BigInteger.valueOf(p);
            if (x.mod(P).equals(BigInteger.ZERO)) {
                out.add(p);
                while (x.mod(P).equals(BigInteger.ZERO)) x = x.divide(P);
            }
        }
        if (x.compareTo(BigInteger.ONE) > 0 && x.bitLength() <= 31) {
            out.add(x.intValue());
        }
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

    /** sqrt(2N) mod q=(∏ p)^2: Tonelli–Shanks (mod p) + Hensel to p^2, CRT-combine. */
    private BigInteger sqrtN2ModCompositeSquare(BigInteger q, ArrayList<Integer> primes) {
        BigInteger r = BigInteger.ZERO;   // current residue
        BigInteger m = BigInteger.ONE;    // current modulus

        for (int p : primes) {
            if (p <= 2) return null;
            BigInteger P = BigInteger.valueOf(p);

            // a ≡ 2N (mod p^2)
            BigInteger p2 = P.multiply(P);
            BigInteger a = N2.mod(p2);

            // sqrt mod p
            long[] roots = MathUtils.ressol(p, a.mod(P).longValue());
            long r0 = -1;
            for (long v : roots) { if (v >= 0) { r0 = v; break; } }
            if (r0 < 0) return null;

            BigInteger s = BigInteger.valueOf(r0);      // s ≡ sqrt(a) (mod p)

            // Hensel lift to p^2:
            // Find k mod p: 2 s k ≡ (a - s^2)/p (mod p)
            BigInteger s2 = s.multiply(s).mod(p2);
            BigInteger diff = a.subtract(s2);
            if (diff.signum() < 0) diff = diff.add(p2);
            BigInteger rhs = diff.divide(P).mod(P);
            BigInteger inv2s = s.shiftLeft(1).mod(P).modInverse(P);
            BigInteger k = rhs.multiply(inv2s).mod(P);
            BigInteger sLift = s.add(k.multiply(P)).mod(p2); // sqrt mod p^2

            // CRT combine with current (r mod m) and (sLift mod p^2)
            BigInteger g = m.gcd(p2);
            if (!g.equals(BigInteger.ONE)) return null;
            BigInteger inv = m.modInverse(p2);
            BigInteger t = sLift.subtract(r).mod(p2);
            BigInteger coeff = t.multiply(inv).mod(p2);
            r = r.add(m.multiply(coeff));
            m = m.multiply(p2);
        }

        return r.mod(m);
    }

    /** Re-center c near floor(sqrt(2N)) to minimize |2N - c^2|; try k, k±1. */
    private static BigInteger recenterC(BigInteger c0, BigInteger Q, BigInteger twoN) {
        BigInteger s = SqrRoot.bigIntSqRootFloor(twoN);
        BigInteger num = s.subtract(c0);
        BigInteger[] dq = num.divideAndRemainder(Q);
        BigInteger k = dq[0];
        if (dq[1].abs().shiftLeft(1).compareTo(Q) >= 0) {
            k = dq[1].signum() >= 0 ? k.add(BigInteger.ONE) : k.subtract(BigInteger.ONE);
        }
        BigInteger bestC = c0.add(k.multiply(Q));
        BigInteger bestGap = bestC.multiply(bestC).subtract(twoN).abs();

        for (int delta = -1; delta <= 1; delta++) {
            BigInteger cand = c0.add(k.add(BigInteger.valueOf(delta)).multiply(Q));
            BigInteger gap = cand.multiply(cand).subtract(twoN).abs();
            if (gap.compareTo(bestGap) < 0) { bestGap = gap; bestC = cand; }
        }
        return bestC;
    }
}
