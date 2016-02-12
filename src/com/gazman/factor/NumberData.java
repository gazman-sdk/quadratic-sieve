package com.gazman.factor;

import java.math.BigInteger;

/**
 * Created by Ilya Gazman on 2/2/2016.
 */
public class NumberData {
    private BigInteger originalX;
    private BigInteger originalY;
    private int version;
    private static int generalVersion;
    private final BigInteger N;
    private BigInteger[] primeBase;
    private BigInteger root;
    private BigInteger y;

    public NumberData(BigInteger N, BigInteger[] primeBase, BigInteger root) {
        this.N = N;
        this.primeBase = primeBase;
        this.root = root;
    }

    public static void upgrade() {
        generalVersion++;
    }

    public int sieve(int primeIndex, long position) {
        if (version < generalVersion) {
            originalX = root.add(BigInteger.valueOf(position));
            originalY = originalX.pow(2).subtract(N);
            y = originalY;
            version = generalVersion;
        }

        int power = 0;
        BigInteger prime = primeBase[primeIndex];
        do {
            if(!y.mod(prime).equals(BigInteger.ZERO)){
                throw new IllegalStateException(y + " is not divided by " + prime);
            }
            y = y.divide(prime);
            power++;
        } while (y.mod(prime).equals(BigInteger.ZERO));

        return y.equals(BigInteger.ONE) ? power : -power;
    }

    public VectorData getData() {
        return new VectorData(originalX, originalY);
    }
}
