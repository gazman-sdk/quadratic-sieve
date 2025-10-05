package com.gazman.factor;

import java.math.BigInteger;
import java.util.*;

/**
 * K-way (up to k=5) large-prime hyper-merge.
 *
 * Stores partial relations keyed by the set (mod 2) of leftover primes (size 1..4).
 * When a new partial arrives, we greedily cancel subsets already seen until we
 * either produce a B-smooth relation (empty set) or we can store the residue.
 *
 * Thread-safe, low-contention (single synchronized add).
 *
 * NOTE:
 * - Each VectorData must carry x and y so multiplying two partials preserves
 *   X^2 ≡ Y (mod N) without knowing the factors of N.
 * - Leftover primes here are NOT added as columns; we only emit B-smooth merges.
 */
public class LargePrimeKMerger {

    private final BigInteger N;

    // Maps by set size: 1..4
    private final Map<String, VectorData> map1 = new HashMap<>();
    private final Map<String, VectorData> map2 = new HashMap<>();
    private final Map<String, VectorData> map3 = new HashMap<>();
    private final Map<String, VectorData> map4 = new HashMap<>();

    public LargePrimeKMerger(BigInteger N) {
        this.N = N;
    }

    /**
     * Attempt to merge a partial with leftover primes 'S' (odd exponents).
     * @param primesOdd list of DISTINCT primes (odd multiplicity) length 1..4
     * @param cur the relation (must have x,y set)
     * @return a B-smooth VectorData if a full cancellation was found; otherwise null and we keep residue.
     */
    public synchronized VectorData add(long[] primesOdd, VectorData cur) {
        if (primesOdd == null || primesOdd.length == 0) {
            return cur; // already B-smooth
        }
        long[] S = primesOdd.clone();
        Arrays.sort(S);

        // Greedy subset-cancel loop. Max 14 subsets for size 4 → cheap.
        while (true) {
            boolean improved = false;

            int n = S.length;
            // Enumerate non-empty proper subsets as bitmasks 1..(2^n-2)
            int maxMask = (1 << n) - 1;
            for (int mask = 1; mask < maxMask; mask++) {
                if (mask == maxMask) continue; // skip full set
                long[] subset = project(S, mask);
                String key = key(subset);
                Map<String, VectorData> bucket = bucket(subset.length);
                VectorData prev = bucket.get(key);
                if (prev != null) {
                    // merge and reduce S ← S Δ subset  (since we track mod 2)
                    cur = merge(cur, prev);
                    S = symmetricDiff(S, subset);
                    improved = true;
                    // Remove the used subset (spent)
                    bucket.remove(key);
                    break; // restart with new S
                }
            }

            if (!improved) {
                // Try full-set collision: same residue set seen → cancel fully
                String fullKey = key(S);
                Map<String, VectorData> bucket = bucket(S.length);
                VectorData prev = bucket.get(fullKey);
                if (prev != null) {
                    bucket.remove(fullKey);
                    return merge(cur, prev); // B-smooth
                }
                // Store residue and return
                bucket.put(fullKey, cur);
                return null;
            }

            // Continue until no subset is found; if S empties, we’re done
            if (S.length == 0) return cur; // B-smooth
        }
    }

    private Map<String, VectorData> bucket(int size) {
        return switch (size) {
            case 1 -> map1;
            case 2 -> map2;
            case 3 -> map3;
            case 4 -> map4;
            default -> throw new IllegalArgumentException("Unsupported residue size " + size);
        };
    }

    private static long[] project(long[] arr, int mask) {
        int n = Integer.bitCount(mask);
        long[] out = new long[n];
        for (int i = 0, j = 0; i < arr.length; i++) {
            if (((mask >> i) & 1) == 1) out[j++] = arr[i];
        }
        return out;
    }

    private static String key(long[] arr) {
        // sorted already
        StringBuilder sb = new StringBuilder(arr.length * 12);
        for (int i = 0; i < arr.length; i++) {
            if (i > 0) sb.append('#');
            sb.append(arr[i]);
        }
        return sb.toString();
    }

    private static long[] symmetricDiff(long[] a, long[] b) {
        // a and b are sorted, unique
        ArrayList<Long> out = new ArrayList<>(a.length + b.length);
        int i = 0, j = 0;
        while (i < a.length || j < b.length) {
            long x = (i < a.length) ? a[i] : Long.MAX_VALUE;
            long y = (j < b.length) ? b[j] : Long.MAX_VALUE;
            if (x == y) { i++; j++; } // cancel (mod 2)
            else if (x < y) { out.add(x); i++; }
            else { out.add(y); j++; }
        }
        long[] r = new long[out.size()];
        for (int k = 0; k < r.length; k++) r[k] = out.get(k);
        return r;
    }

    private VectorData merge(VectorData a, VectorData b) {
        // XOR parities
        BitSetClone bits = new BitSetClone(a.vector);
        bits.xor(b.vector);

        VectorData m = new VectorData(bits.bitSet, -1);
        // Multiply x and y. x is always reduced mod N in your solver path later,
        // but storing reduced now keeps numbers smaller.
        if (a.x != null && b.x != null) {
            m.x = a.x.multiply(b.x).mod(N);
        }
        if (a.y != null && b.y != null) {
            m.y = a.y.multiply(b.y);
        }
        return m;
    }

    // Tiny helper: clone to avoid mutating original BitSets
    private static class BitSetClone {
        final java.util.BitSet bitSet;
        BitSetClone(java.util.BitSet src) { this.bitSet = (java.util.BitSet) src.clone(); }
        void xor(java.util.BitSet other) { this.bitSet.xor(other); }
    }
}
