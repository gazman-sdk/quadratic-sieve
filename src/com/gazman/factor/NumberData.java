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
    private int divisions;

    public NumberData(BigInteger N, BigInteger[] primeBase, BigInteger root) {
        this.N = N;
        this.primeBase = primeBase;
        this.root = root;
    }

    public static void upgrade() {
        generalVersion++;
    }

    public boolean sieve(int primeIndex, long position) {
        if (version < generalVersion) {
            divisions = 0;
            originalX = root.add(BigInteger.valueOf(position));
            originalY = originalX.pow(2).subtract(N);
            y = originalY;
            version = generalVersion;
        }

        BigInteger prime = primeBase[primeIndex];
        do {
            y = y.divide(prime);
            divisions++;
        } while (y.mod(prime).equals(BigInteger.ZERO));

        return y.equals(BigInteger.ONE);
    }

    public BigInteger getRemindY(){
        return version < generalVersion || y.equals(BigInteger.ONE)? null : y;
    }

    public VectorData getData() {
        return new VectorData(originalX, originalY);
    }

    public int getDivisions() {
        return divisions;
    }
}
