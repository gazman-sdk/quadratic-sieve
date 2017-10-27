package com.gazman.factor;

import java.lang.reflect.Array;
import java.math.BigInteger;
import java.util.List;

/**
 * Created by Ilya Gazman on 1/11/2016.
 */
public class Logger {

    private static final long startingTime = System.currentTimeMillis();
    final BigInteger one = BigInteger.ONE;

    private final boolean logsEnabled = true;


    private final String emptyChar = (new char[1])[0] + "";

    protected synchronized void log(Object... params) {
        if (!logsEnabled) {
            return;
        }
        if (params.length > 0) {
            logTime();
            logInLine(params);
        }
        System.out.println();
    }

    private void logTime() {
        long time = System.currentTimeMillis() - startingTime;
        long seconds = time / 1000;
        long milliseconds = time - seconds * 1000;
        String millisecondsString = milliseconds + "";
        switch (millisecondsString.length()) {
            case 1:
                millisecondsString = "00" + millisecondsString;
                break;
            case 2:
                millisecondsString = "0" + millisecondsString;
                break;
        }
        System.out.print(seconds + "." + millisecondsString + ": ");
    }

    private void logInLine(Object... params) {
        if (!logsEnabled) {
            return;
        }
        for (Object param : params) {
            if (param instanceof Number) {
                printNumber(param);
            } else if (param instanceof Boolean) {
                System.out.print((Boolean) param ? 1 : 0);
            } else if (param.getClass().isArray()) {
                System.out.print("[");
                for (int i = 0; i < Array.getLength(param); i++) {
                    System.out.print(Array.get(param, i) + " ");
                }
                System.out.print("]");
            } else if (param instanceof List) {
                System.out.print("[");
                List list = (List) param;
                for (int i = 0; i < list.size(); i++) {
                    Object o = list.get(i);
                    System.out.print(o + (i < list.size() - 1 ? ", " : ""));
                }
                System.out.print("]");
            } else {
                System.out.print(param + " ");
            }
        }
    }

    private void printNumber(Object value) {
        int minLength = 2;
        String out = value.toString() + " ";
        out = updateLength(out, minLength);
        System.out.print(out);
    }

    private String updateLength(String input, int minLength) {
        if (input.length() < minLength) {
            input = input + new String(new char[minLength - input.length()]).replaceAll(emptyChar, " ");
        }
        return input;
    }
}
