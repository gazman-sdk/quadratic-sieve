package com.gazman.factor.matrix;

import com.gazman.factor.BigPrimesList;
import com.gazman.factor.VectorData;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.LinkedList;

/**
 * Single-large-prime merging for k=5.
 * Produces B-smooth relations by merging two 1-LP partials sharing the same big prime.
 * For each merged relation:
 *   x = Π A(X_i) mod N,  A(X)=f^2 s^2 (sX - t)(sX + t)  (mod N), with X = M + x_i
 *   y = Π Q5(X_i)  (integer)
 * and parity vector = XOR.
 */
public class VectorsShrinker {

    private BigInteger N;
    private int bigPrimesIndex;

    // k=5 polynomial + center
    private BigInteger f, s, t, c, MB;

    public void initK5(int biggestPrimeIndex, BigInteger N,
                       BigInteger f, BigInteger s, BigInteger t, BigInteger c,
                       BigInteger MB) {
        this.bigPrimesIndex = biggestPrimeIndex;
        this.N = N;
        this.f = f; this.s = s; this.t = t; this.c = c; this.MB = MB;
    }

    @SuppressWarnings("unchecked")
    public ArrayList<VectorData> shrink(ArrayList<VectorData> bSmoothVectors, BigPrimesList bigPrimesList) {
        bSmoothVectors = (ArrayList<VectorData>) bSmoothVectors.clone();
        ArrayList<LinkedList<VectorData>> bigPrimes = bigPrimesList.getBigPrimes();

        // Exact pairs first
        for (int i = bigPrimes.size() - 1; i >= 0; i--) {
            LinkedList<VectorData> vectorDatas = bigPrimes.get(i);
            if (vectorDatas.size() == 2) {
                bigPrimes.remove(i);
                VectorData a = vectorDatas.get(0);
                VectorData b = vectorDatas.get(1);
                mergeInto(a, b);
                bSmoothVectors.add(a);
            }
        }

        // Residual big primes: assign synthetic indices
        for (LinkedList<VectorData> bigPrimeList : bigPrimes) {
            boolean updateIndex = false;
            boolean firstVector = true;
            int bigPrimeIndex = -1;

            for (VectorData vd : bigPrimeList) {
                if (firstVector && vd.bigPrimeIndex == -1) {
                    firstVector = false;
                    updateIndex = true;
                }
                if (vd.bigPrimeIndex == -1) {
                    if (bigPrimeIndex == -1) bigPrimeIndex = this.bigPrimesIndex;
                    vd.bigPrimeIndex = bigPrimeIndex;
                    vd.vector.set(bigPrimeIndex);
                } else {
                    bigPrimeIndex = vd.bigPrimeIndex;
                }
                bSmoothVectors.add(vd);
            }
            if (updateIndex) this.bigPrimesIndex++;
        }

        return bSmoothVectors;
    }

    private void mergeInto(VectorData result, VectorData other) {
        result.vector.xor(other.vector);

        BigInteger X1 = MB.add(BigInteger.valueOf(result.position));
        BigInteger X2 = MB.add(BigInteger.valueOf(other.position));

        BigInteger A1 = A_of(X1);
        BigInteger A2 = A_of(X2);
        BigInteger Q1 = Q5_of(X1);
        BigInteger Q2 = Q5_of(X2);

        result.x = A1.multiply(A2).mod(N);
        result.y = Q1.multiply(Q2);
    }

    // A(X) = f^2 s^2 (sX - t)(sX + t)  (mod N)
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
