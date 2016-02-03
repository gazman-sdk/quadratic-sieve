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
        long exponent = p.subtract(BigInteger.ONE).divide(tow).longValue();
        return modularExponentiation(x.longValue(), exponent, p.longValue()) == 1;
    }

    // based on http://math.stackexchange.com/a/453108/101178
    public static long modularExponentiation(long value, long exponent, long mod) {
        long result = 1;
        while (exponent > 0) {
            if ((exponent & 1) == 1) {
                value = value % mod;
                result = (result * value) % mod;
                result = result % mod;
            }
            exponent = exponent >> 1;
            value = value % mod;
            value = (value * value) % mod;
            value = value % mod;
        }
        return result;
    }

}
