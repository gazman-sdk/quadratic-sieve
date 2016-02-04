package com.gazman.factor;

import java.math.BigInteger;

/**
 * Created by Ilya Gazman on 2/2/2016.
 */
public class VectorNumber {
    private BigInteger originalX;
    private BigInteger originalY;
    private int version;
    private static int generalVersion;
    private final BigInteger N;
    private BigInteger[] primeBase;
    private BigInteger root;
    private BigInteger y;

    public VectorNumber(BigInteger N, BigInteger[] primeBase, BigInteger root) {
        this.N = N;
        this.primeBase = primeBase;
        this.root = root;
    }

    public static void upgrade() {
        generalVersion++;
    }

    public boolean sieve(int primeIndex, long position) {
        if (version < generalVersion) {
            originalX = root.add(BigInteger.valueOf(position)).pow(2);
            originalY = originalX.subtract(N);
            y = originalY;
            version = generalVersion;
        }

        BigInteger prime = primeBase[primeIndex];
        do {
            if(!y.mod(prime).equals(BigInteger.ZERO)){
                break;
            }
            y = y.divide(prime);
        } while (y.mod(prime).equals(BigInteger.ZERO));

        return y.equals(BigInteger.ONE);
    }

    public VectorData getData() {
        return new VectorData(originalX, originalY);
    }
}
