package com.gazman;

import com.gazman.factor.Logger;
import com.gazman.search.FSTSearch;

import java.io.FileDescriptor;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.math.BigInteger;
import java.nio.charset.StandardCharsets;
import java.util.Random;

public class SearchMain extends Logger {
    private static final Random random = new Random(12);

    static void main(String[] args) {
        forceUtf8Stdout();
        new SearchMain().init();
    }

    public static void forceUtf8Stdout()  {
        System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out), true, StandardCharsets.UTF_8));
        System.setErr(new PrintStream(new FileOutputStream(FileDescriptor.err), true, StandardCharsets.UTF_8));
    }


    private void init() {
        int length = 90;
        BigInteger a = BigInteger.probablePrime(length + 1, random);
        BigInteger b = BigInteger.probablePrime(length - 1, random);

        BigInteger input = a.multiply(b);
        log("Input N:", input);
        log("N bits:", input.bitLength());
        log("---------");
        log();

        // ------ Try the robust k=5 search ------
        try {
            // A rough estimate for the sieve interval size for scoring purposes
            long sieveM = 1L << (input.bitLength() / 5);
            FSTSearch fst = new FSTSearch(input);
            FSTSearch.Poly best = fst.findBestPolynomial(sieveM);

            log();
            log("------ k=5 Search Successful ------");
            log("Best polynomial found with score:", String.format("%.2f", best.score));
            log("  f =", best.f);
            log("  s =", best.s, "(" + best.s.bitLength() + " bits)");
            log("  t =", best.t, "(" + best.t.bitLength() + " bits)");
            log("  c =", best.c, "(" + best.c.bitLength() + " bits)");

            // Verify the predicted bit sizes at the edge of the sieve interval
            int targetBits = input.bitLength() / 5;
            BigInteger M_BI = BigInteger.valueOf(sieveM);
            BigInteger fs2 = best.f.multiply(best.s).multiply(best.s);
            int l1_bits = best.s.multiply(M_BI).subtract(best.t).abs().bitLength();
            int l2_bits = best.s.multiply(M_BI).add(best.t).abs().bitLength();
            int l3_bits = fs2.multiply(M_BI).subtract(best.c).abs().bitLength();
            int l4_bits = fs2.multiply(M_BI).add(best.c).abs().bitLength();
            int const_bits = best.f.multiply(best.s).bitLength() * 2;

            log();
            log("Predicted max bit sizes at sieve edge (M=" + sieveM + "):");
            log("  Target:    ", targetBits);
            log("  L1(sx-t):  ", l1_bits);
            log("  L2(sx+t):  ", l2_bits);
            log("  L3(fs²x-c):", l3_bits);
            log("  L4(fs²x+c):", l4_bits);
            log("  Const(f²s²):", const_bits);

        } catch (Exception ex) {
            log();
            log("------ k=5 Search FAILED ------");
            log(ex.getMessage());
        }
    }
}