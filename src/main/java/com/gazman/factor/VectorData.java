package com.gazman.factor;

import java.math.BigInteger;
import java.util.BitSet;

/**
 * Created by Ilya Gazman on 1/30/2016.
 */
public class VectorData {

    public final BitSet vector;
    public final long position;
    public int bigPrimeIndex = -1;

    public BigInteger x;
    public BigInteger y;

    VectorData(BitSet vector, long position) {
        this.vector = vector;
        this.position = position;
    }

    @Override
    public String toString() {
        return vector.toString();
    }

    @Override
    public int hashCode() {
        return vector.toString().hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        return obj instanceof VectorData && vector.toString().equals(((VectorData) obj).vector.toString());
    }
}
