package com.gazman.math;

import java.math.BigInteger;

/**
 * Source: http://stackoverflow.com/questions/4407839/how-can-i-find-the-square-root-of-a-java-biginteger
 */

@SuppressWarnings("SpellCheckingInspection")
public class SqrRoot {
    private static final BigInteger TWO = BigInteger.valueOf(2L);

    public static BigInteger bigIntSqRootFloor(BigInteger value)
            throws IllegalArgumentException {
        if (checkTrivial(value)) {
            return value;
        }
        if (value.bitLength() < 64) { // Can be cast to long
            double sqrt = Math.sqrt(value.longValue());
            return BigInteger.valueOf((long) sqrt);
        }
        // starting with y = value / 2 avoids magnitude issues with value squared
        BigInteger y = TWO.pow(value.bitLength() / 2);
        BigInteger result = value.divide(y);
        while (!y.subtract(result).abs().equals(BigInteger.ONE) && !y.subtract(result).equals(BigInteger.ZERO)) {
            y = result.add(y).divide(TWO);
            result = value.divide(y);
        }
        return y;
    }

    public static BigInteger bigIntSqRootCeil(BigInteger x)
            throws IllegalArgumentException {
        BigInteger y = bigIntSqRootFloor(x);
        if (x.compareTo(y.multiply(y)) == 0) {
            return y;
        }
        return y.add(BigInteger.ONE);
    }

    private static boolean checkTrivial(BigInteger x) {
        if (x == null) {
            throw new NullPointerException("x can't be null");
        }
        if (x.compareTo(BigInteger.ZERO) < 0) {
            throw new IllegalArgumentException("Negative argument.");
        }

        return x.equals(BigInteger.ZERO) || x.equals(BigInteger.ONE);
    }
}
