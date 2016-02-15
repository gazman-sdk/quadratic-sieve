package com.gazman.factor;

import java.util.HashMap;
import java.util.LinkedList;

/**
 * Created by Ilya Gazman on 2/15/2016.
 */
public class BigPrimesList {

    private static final Double PRECISION = Math.pow(10, 8);
    private LinkedList<LinkedList<VectorData>> primes = new LinkedList<>();
    private HashMap<Long, LinkedList<VectorData>> primesHash = new HashMap<>();
    private int primesFound;

    public void add(Double primeLog, VectorData bigPrime){
        long prime = Math.round(Math.pow(Math.E, primeLog));

        LinkedList<VectorData> vectorDatas = primesHash.get(prime);
        if(vectorDatas == null){
            vectorDatas = new LinkedList<>();
            primesHash.put(prime, vectorDatas);
            primes.add(vectorDatas);
        }
        vectorDatas.add(bigPrime);
        if(vectorDatas.size() > 1){
            primesFound++;
        }
    }

    public int getPrimesFound() {
        return primesFound;
    }

    public LinkedList<LinkedList<VectorData>> getBigPrimes(){
        LinkedList<LinkedList<VectorData>> primeList = new LinkedList<>();
        for (LinkedList<VectorData> prime : primes) {
            if(prime.size() > 1){
                primeList.add(prime);
            }
        }

        return primeList;
    }
}
