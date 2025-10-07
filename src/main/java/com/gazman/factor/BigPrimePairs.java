package com.gazman.factor;

import java.math.BigInteger;
import java.util.BitSet;
import java.util.HashMap;

/**
 * Collects 2-large-prime partials keyed by an unordered pair (p,q).
 * For k=5, merged VectorData carries:
 *   x = Π A(X_i) mod N,     A(X)=f^2 s^2 (sX - t)(sX + t)  (mod N), with X = M + x_i
 *   y = Π Q5(X_i) as integers,
 * and parity = XOR.
 */
public class BigPrimePairs {

    private final HashMap<String, VectorData> pairs = new HashMap<>();
    private final BigInteger N;
    private final BigInteger f, s, t, c, MB;

    public BigPrimePairs(BigInteger N, BigInteger f, BigInteger s, BigInteger t, BigInteger c, BigInteger MB) {
        this.N = N;
        this.f = f; this.s = s; this.t = t; this.c = c; this.MB = MB;
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
        // Merge
        BitSet vec = (BitSet) prev.vector.clone();
        vec.xor(vd.vector);
        VectorData merged = new VectorData(vec, -1);

        BigInteger X1 = MB.add(BigInteger.valueOf(prev.position));
        BigInteger X2 = MB.add(BigInteger.valueOf(vd.position));

        BigInteger A1 = A_of(X1);
        BigInteger A2 = A_of(X2);
        BigInteger Q1 = Q5_of(X1);
        BigInteger Q2 = Q5_of(X2);

        merged.x = A1.multiply(A2).mod(N);
        merged.y = Q1.multiply(Q2);
        return merged;
    }

    // A(X) = f^2 s^2 (sX - t)(sX + t) (mod N)
    private BigInteger A_of(BigInteger X) {
        BigInteger sX = s.multiply(X);
        BigInteger term = sX.subtract(t).multiply(sX.add(t));
        return f.multiply(f).multiply(s).multiply(s).multiply(term).mod(N);
    }

    // Q5(X) = (f^2 s^2) * moving part
    private BigInteger Q5_of(BigInteger X) {
        BigInteger sX   = s.multiply(X);
        BigInteger fs2X = f.multiply(s).multiply(s).multiply(X);
        BigInteger mov  = sX.subtract(t).multiply(sX.add(t))
                .multiply(fs2X.subtract(c))
                .multiply(fs2X.add(c));
        return f.multiply(f).multiply(s).multiply(s).multiply(mov);
    }
}
