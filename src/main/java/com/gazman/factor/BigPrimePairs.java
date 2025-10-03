package com.gazman.factor;

import java.math.BigInteger;
import java.util.BitSet;
import java.util.HashMap;

/**
 * Collects 2-large-prime partials keyed by an unordered pair (p,q).
 * When a second relation with the same (p,q) arrives, returns
 * a merged B-smooth VectorData (parity = XOR, x/y multiplied).
 */
public class BigPrimePairs {

    private final HashMap<String, VectorData> pairs = new HashMap<>();
    private final BigInteger N;
    private final BigInteger root;

    public BigPrimePairs(BigInteger root, BigInteger N) {
        this.root = root;
        this.N = N;
    }

    private static String key(long a, long b) {
        long x = Math.min(a, b), y = Math.max(a, b);
        return x + "#" + y;
    }

    public synchronized VectorData add(long p, long q, VectorData vd) {
        String k = key(p, q);
        VectorData prev = pairs.remove(k);
        if (prev == null) {
            pairs.put(k, vd);
            return null;
        }
        // Merge to a B-smooth relation
        BitSet vec = (BitSet) prev.vector.clone();
        vec.xor(vd.vector);
        VectorData merged = new VectorData(vec, -1);

        BigInteger x1 = root.add(BigInteger.valueOf(prev.position));
        BigInteger x2 = root.add(BigInteger.valueOf(vd.position));
        BigInteger y1 = x1.multiply(x1).subtract(N);
        BigInteger y2 = x2.multiply(x2).subtract(N);

        merged.x = x1.multiply(x2);
        merged.y = y1.multiply(y2);
        return merged;
    }
}
