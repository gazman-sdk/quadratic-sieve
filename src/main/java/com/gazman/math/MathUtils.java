package com.gazman.math;

import java.math.BigInteger;

/**
 * Created by Ilya Gazman on 2/1/2016.
 */
public class MathUtils {
    public static boolean isRootInQuadraticResidues(BigInteger n, BigInteger p) {
        BigInteger tow = BigInteger.valueOf(2);
        BigInteger x = n.mod(p);
        if (p.equals(tow)) {
            return x.mod(tow).equals(BigInteger.ONE);
        }
        BigInteger exponent = p.subtract(BigInteger.ONE).divide(tow);
        return x.modPow(exponent, p).equals(BigInteger.ONE);
    }

    public static double log(double x, double base) {
        return Math.log(x) / Math.log(base);
    }

    /**
     * Tonelli–Shanks algorithm implemented by Stefan Buettcher
     *
     * @param prime the prime modular
     * @param n     quadratic residue to find solutions for
     * @return array of longs with size 2, if there is no solution the value will be -1
     * @see "http://www.stefan.buettcher.org/cs/factorization/index.html"
     */
    public static long[] ressol(long prime, long n) {
        long k, x;
        BigInteger bigN = BigInteger.valueOf(n);
        BigInteger bigPrime = BigInteger.valueOf(prime);

        long result[] = new long[2];
        result[0] = -1;
        result[1] = -1;

        if (prime == 2) {
            result[0] = n % 2;
            result[1] = -1;
            return result;
        }

        if (prime % 4 == 3) {
            k = (prime / 4);

            x = modPowLong(bigN, k + 1, bigPrime) % prime;
            result[0] = x;
            result[1] = (prime - x);
            return result;
        }

        if (prime % 8 == 5) {
            k = (prime / 8);
            x = modPowLong(bigN, 2 * k + 1, bigPrime);
            if (x == 1) {
                x = modPowLong(bigN, k + 1, bigPrime);
                result[0] = x;
                result[1] = (prime - x);
                return result;
            }
            if (x == prime - 1) {
                x = modPowLong(BigInteger.valueOf(4 * n), k + 1, bigPrime);
                x = (x * (prime + 1) / 2) % prime;
                result[0] = x;
                result[1] = (prime - x);
                return result;
            }
        }

        long h = 13;
        do {
            h += 2;
        }
        while (isRootInQuadraticResidues(BigInteger.valueOf(h * h - 4 * n),
                bigPrime));

        k = (prime + 1) / 2;
        x = v_(k, h, n, prime);
        if (x < 0) {
            x += prime;
        }
        x = (x * k) % prime;
        result[0] = x;
        result[1] = (prime - x);

        return result;
    }

    private static long modPowLong(BigInteger n, long exponent, BigInteger modulo) {
        return n.modPow(BigInteger.valueOf(exponent), modulo).longValue();
    }

    private static long v_(long j, long h, long n, long p) {
        long b[] = new long[64];
        long m = n;
        long v = h;
        long w = (h * h - 2 * m) % p;
        long x;
        int k, t;
        t = 0;
        while (j > 0) {
            b[++t] = j % 2;
            j /= 2;
        }
        for (k = t - 1; k >= 1; k--) {
            x = (v * w - h * m) % p;
            v = (v * v - 2 * m) % p;
            w = (w * w - 2 * n * m) % p;
            m = m * m % p;
            if (b[k] == 0)
                w = x;
            else {
                v = x;
                m = n * m % p;
            }
        }
        return v;
    }

}
