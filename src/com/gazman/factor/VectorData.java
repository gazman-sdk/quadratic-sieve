package com.gazman.factor;

import java.math.BigInteger;

/**
 * Created by Ilya Gazman on 1/30/2016.
 */
public class VectorData {
    public BigInteger x,y;
    public boolean[] vector;

    public VectorData(BigInteger x, BigInteger y) {
        this.x = x;
        this.y = y;
    }

    public void buildVector(BigInteger primeBase[]) {
        if(vector != null){
            return;
        }
        BigInteger y = this.y;
        vector = new boolean[primeBase.length];
        for (int i = 0; i < primeBase.length && !y.equals(BigInteger.ONE); i++) {
            BigInteger prime = primeBase[i];
            while (y.mod(prime).equals(BigInteger.ZERO)){
                y = y.divide(prime);
                vector[i] = !vector[i];
            }
        }
    }
}
