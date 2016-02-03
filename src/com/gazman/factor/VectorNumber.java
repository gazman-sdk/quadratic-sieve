package com.gazman.factor;

import java.math.BigInteger;

/**
 * Created by Ilya Gazman on 2/2/2016.
 */
public class VectorNumber {
    private BigInteger originalX;
    private int version;
    private static int generalVersion;
    private final BigInteger N;
    private BigInteger[] primeBase;
    private BigInteger root;
    private VersionBoolean[] vector;
    private BigInteger originalY;
    private BigInteger y;

    public VectorNumber(BigInteger N, BigInteger[] primeBase, BigInteger root) {
        this.N = N;
        this.primeBase = primeBase;
        this.root = root;
        vector = new VersionBoolean[primeBase.length];
        for (int i = 0; i < vector.length; i++) {
            vector[i] = new VersionBoolean();
        }
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
            vector[primeIndex].not();
        } while (y.mod(prime).equals(BigInteger.ZERO));

        return y.equals(BigInteger.ONE);
    }

    public VectorData buildVector(){
        VectorData vectorData = new VectorData();
        vectorData.x = originalX;
        vectorData.y = originalY;
        vectorData.vector = new boolean[this.vector.length];
        for (int i = 0; i < vectorData.vector.length; i++) {
            vectorData.vector[i] = this.vector[i].getValue();
        }

        return vectorData;
    }

    private class VersionBoolean {
        private int version;
        private boolean value;

        public boolean getValue() {
            return version >= generalVersion && value;
        }

        public void not() {
            value = !getValue();
            version = generalVersion;
        }

        public void setValue(boolean value) {
            version = generalVersion;
            this.value = value;
        }
    }
}
