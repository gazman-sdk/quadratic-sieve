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

        // Run standalone search (self-logging)
        ABSearch.setStaticN(input);
        ABSearch search = new ABSearch(input);

        int msBudget = 20_000;     // 20 seconds
        int topK     = 6;          // show top 6
        int pmax     = 2_000_000;  // prime ceiling for pool
        Integer xMax = null;       // let search choose default

        List<ABSearch.Poly> result = search.runFixedTime(msBudget, topK, pmax, xMax);
        // ABSearch prints a summary itself; nothing else required here.
    }
}
