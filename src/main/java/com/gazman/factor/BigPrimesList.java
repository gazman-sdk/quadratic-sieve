package com.gazman.factor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;

/**
 * Created by Ilya Gazman on 2/15/2016.
 */
public class BigPrimesList {

    private final LinkedList<LinkedList<VectorData>> primes = new LinkedList<>();
    private final HashMap<Long, LinkedList<VectorData>> primesHash = new HashMap<>();
    private int primesFound;

    void add(long prime, VectorData bigPrime) {
        LinkedList<VectorData> vectorDatas = primesHash.get(prime);
        if (vectorDatas == null) {
            vectorDatas = new LinkedList<>();
            primesHash.put(prime, vectorDatas);
            primes.add(vectorDatas);
        }
        vectorDatas.add(bigPrime);
        if (vectorDatas.size() > 1) {
            primesFound++;
        }
    }

    int getPrimesFound() {
        return primesFound;
    }

    public ArrayList<LinkedList<VectorData>> getBigPrimes() {
        ArrayList<LinkedList<VectorData>> primeList = new ArrayList<>();
        for (LinkedList<VectorData> prime : primes) {
            if (prime.size() > 1) {
                primeList.add(prime);
            }
        }

        return primeList;
    }
}
