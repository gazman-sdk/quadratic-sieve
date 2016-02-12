package com.gazman;

import com.gazman.factor.*;
import com.gazman.math.MathUtils;
import com.gazman.math.SqrRoot;

import java.math.BigDecimal;
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
//        a = BigInteger.valueOf(139);
//        b = BigInteger.valueOf(191);
//        a = BigInteger.valueOf(900289138583L);
//        b = BigInteger.valueOf(633491976539L);

//        a = BigInteger.valueOf(648336205039L);
//        b = BigInteger.valueOf(985468814843L);


        BigInteger input = a.multiply(b);
        log(a, b);
        log(input, input.toString().length());
        log("---------");
        log();
        new QuadraticThieve().factor(input);
    }
}
