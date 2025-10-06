package com.gazman;

import com.gazman.factor.Logger;
import com.gazman.search.ABSearch;

import java.math.BigInteger;
import java.util.List;
import java.util.Random;

public class SearchMain extends Logger {

    private static final Random random = new Random(12);

    static void main() {
        new SearchMain().init();
    }

    private void init() {
        int length = 90;
        BigInteger a = BigInteger.probablePrime(length + 1, random);
        BigInteger b = BigInteger.probablePrime(length - 1, random);

        BigInteger input = a.multiply(b);
        log(a, b);
        log(input, input.toString().length());
        log("---------");
        log();

        ABSearch search = new ABSearch(input);

        int msBudget = 20_000;     // 20 seconds
        int topK = 16;         // keep more top results (helps surface richer A)
        int pmax = 5_000_000;  // richer pool → more factors possible in A
        Integer xMax = null;       // let search choose default

        int threads = Math.max(1, Runtime.getRuntime().availableProcessors());
        List<ABSearch.Poly> result = search.runFixedTimeParallel(msBudget, topK, pmax, xMax, threads);
        // Results + summary are logged by ABSearch.
    }
}
