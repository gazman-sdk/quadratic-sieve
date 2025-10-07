package com.gazman;

import com.gazman.factor.Logger;
import com.gazman.factor.QuadraticThieve;

import java.math.BigInteger;
import java.util.Random;

import static com.gazman.SearchMain.forceUtf8Stdout;

public class Main extends Logger {

    private static final Random random = new Random(12);

    static void main() {
        forceUtf8Stdout();
        new Main().init();
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
        new QuadraticThieve(input).start();
    }
}
