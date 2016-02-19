package com.gazman;

import com.gazman.factor.Logger;
import com.gazman.factor.QuadraticThieve;

import java.math.BigInteger;
import java.util.Random;

public class Main extends Logger {

    private static Random random = new Random();

    public static void main(String[] args) {
        new Main().init();
    }

    private void init() {
        int length = 50;
        BigInteger a = BigInteger.probablePrime(length + 1, random);
        BigInteger b = BigInteger.probablePrime(length - 1, random);

//        a = new BigInteger("1209557");
//        b = new BigInteger("280837");

        BigInteger input = a.multiply(b);
//        input = RSA220;
        log(a, b);
        log(input, input.toString().length());
        log("---------");
        log();
        new QuadraticThieve().factor(input);
    }
}
