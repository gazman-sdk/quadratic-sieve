package com.gazman;

import com.gazman.factor.*;

import java.math.BigInteger;
import java.util.Random;

public class Main extends BaseFactor {

    private static Random random = new Random();

    public static void main(String[] args) {
        new Main().init();
    }

    private void init() {
        int length = 20;
        BigInteger a = BigInteger.probablePrime(length, random);
        BigInteger b = BigInteger.probablePrime(length, random);

        BigInteger input = a.multiply(b);
        log(a, b);
        log(input);
        log("---------");
        log();
        new QuadraticThieve().factor(input);
    }
}
