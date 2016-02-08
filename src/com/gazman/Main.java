package com.gazman;

import com.gazman.factor.*;

import java.math.BigInteger;
import java.util.Random;

public class Main extends Logger {

    private static Random random = new Random();

    public static void main(String[] args) {
        new Main().init();
    }

    private void init() {
        int length = 40;
        BigInteger a = BigInteger.probablePrime(length, random);
        BigInteger b = BigInteger.probablePrime(length, random);

//        a = BigInteger.valueOf(179);
//        b = BigInteger.valueOf(149);

        BigInteger input = a.multiply(b);
        log(a, b);
        log(input, input.toString().length());
        log("---------");
        log();
        new QuadraticThieve().factor(input);
    }
}
