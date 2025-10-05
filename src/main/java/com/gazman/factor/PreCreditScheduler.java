package com.gazman.factor;

import com.gazman.math.MathUtils;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Density-aware pre-credit scheduler (CRT over small primes).
 *
 * We stop adding primes when the expected eligible count per block
 * would drop below 'minExpectedHitsPerBlock'. This prevents the
 * “near-zero density” situation.
 */
public class PreCreditScheduler {

    public static final class Summary {
        public final int usedPrimes;              // number of primes in CRT
        public final int modulus;                 // P (fits in int)
        public final int classesCount;            // |R|
        public final double creditBits;           // log2(P)
        public final double density;              // |R| / P
        public final double expectedHitsPerBlock; // density * blockLen

        Summary(int usedPrimes, int modulus, int classesCount, double creditBits, double density, double expectedHitsPerBlock) {
            this.usedPrimes = usedPrimes;
            this.modulus = modulus;
            this.classesCount = classesCount;
            this.creditBits = creditBits;
            this.density = density;
            this.expectedHitsPerBlock = expectedHitsPerBlock;
        }
    }

    private final int[] residuesModP; // sorted residues (t mod P) that pass
    private final int P;              // combined modulus
    private final Summary summary;

    public PreCreditScheduler(BigInteger N,
                              BigInteger root,
                              BigInteger[] primeBase,
                              int[] primeBaseInt,
                              int blockLen,                // sieveVectorBound
                              int targetCreditBits,        // aim, but density floor wins
                              int maxClasses,              // cap |R|
                              int maxModulusBits,          // P <= 2^maxModulusBits
                              int minExpectedHitsPerBlock) // density floor * blockLen
    {
        List<Integer> classes = new ArrayList<>();
        classes.add(0); // 0 mod 1
        int mod = 1;
        int used = 0;
        double credit = 0.0;

        int maxModulus = (maxModulusBits >= 30) ? (1 << 30) : (1 << maxModulusBits);

        for (int i = 0; i < primeBaseInt.length; i++) {
            int p = primeBaseInt[i];
            if (p <= 1) continue;

            long predictedMod = (long) mod * (long) p;
            if (predictedMod > maxModulus) break;

            long[] s = MathUtils.ressol(p, N.mod(BigInteger.valueOf(p)).longValue());
            int sCount = 0;
            for (long v : s) if (v >= 0) sCount++;
            if (sCount == 0) continue; // should not happen for FB primes

            int rootModP = root.mod(BigInteger.valueOf(p)).intValue();
            int[] tResidues = new int[sCount];
            int idx = 0;
            for (long v : s) {
                if (v < 0) continue;
                int r = (int) (v - rootModP);
                r %= p;
                if (r < 0) r += p;
                tResidues[idx++] = r;
            }
            if (idx == 2 && tResidues[0] == tResidues[1]) idx = 1;
            tResidues = Arrays.copyOf(tResidues, idx);

            List<Integer> next = new ArrayList<>(Math.min(maxClasses, classes.size() * tResidues.length));
            int inv = invMod(mod % p, p);
            if (inv == 0) continue;

            for (int a : classes) {
                for (int r : tResidues) {
                    int rhs = r - (a % p);
                    rhs %= p;
                    if (rhs < 0) rhs += p;
                    int k = (int) (((long) rhs * (long) inv) % p);
                    int x = a + mod * k;
                    x %= (mod * p);
                    if (x < 0) x += (mod * p);
                    next.add(x);
                    if (next.size() >= maxClasses) break;
                }
                if (next.size() >= maxClasses) break;
            }

            int newMod = mod * p;
            int newClasses = next.size();
            if (newClasses == 0) break;

            double newDensity = (double) newClasses / (double) newMod;
            double newExpected = newDensity * (double) blockLen;
            if (newExpected < (double) minExpectedHitsPerBlock) {
                break; // stop before overshooting the density floor
            }

            mod = newMod;
            classes = next;
            used++;
            credit = Math.log(mod) / Math.log(2.0);
            if (credit >= targetCreditBits) break;
        }

        P = mod;
        residuesModP = classes.stream().mapToInt(Integer::intValue).sorted().toArray();

        double density = (P == 0) ? 1.0 : ((double) residuesModP.length) / (double) P;
        double expected = density * (double) blockLen;
        summary = new Summary(used, P, residuesModP.length, credit, density, expected);
    }

    public Summary getSummary() { return summary; }

    /** All indices i in [0,len) with (baseT + i) ≡ residue (mod P). */
    public int[] eligibleIndicesForBlock(long baseT, int len) {
        if (P <= 1 || residuesModP.length == 0 || len <= 0) {
            int[] all = new int[len];
            for (int i = 0; i < len; i++) all[i] = i;
            return all;
        }
        int start = (int) Math.floorMod(baseT, P);
        int cap = Math.max(residuesModP.length, (int) Math.ceil(len * summary.density) + residuesModP.length + 8);
        int[] tmp = new int[Math.max(16, cap)];
        int ptr = 0;

        for (int r : residuesModP) {
            int first = r - start;
            if (first < 0) first += P;
            for (int idx = first; idx < len; idx += P) {
                if (ptr == tmp.length) tmp = Arrays.copyOf(tmp, (int) (tmp.length * 1.5) + 8);
                tmp[ptr++] = idx;
            }
        }
        return Arrays.copyOf(tmp, ptr);
    }

    private static int invMod(int a, int p) {
        a %= p;
        if (a < 0) a += p;
        int t = 0, newT = 1;
        int r = p, newR = a;
        while (newR != 0) {
            int q = r / newR;
            int tmpT = t - q * newT; t = newT; newT = tmpT;
            int tmpR = r - q * newR; r = newR; newR = tmpR;
        }
        if (r > 1) return 0;
        if (t < 0) t += p;
        return t;
    }
}
