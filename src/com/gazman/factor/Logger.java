package com.gazman.factor;

import java.lang.reflect.Array;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.List;

/**
 * Created by Ilya Gazman on 1/11/2016.
 */
public class Logger {

    private static long startingTime = System.currentTimeMillis();
    protected BigInteger zero = BigInteger.ZERO;
    protected BigInteger one = BigInteger.ONE;
    protected BigInteger minus_one = BigInteger.valueOf(-1);
    protected BigInteger two = BigInteger.valueOf(2);
    protected BigDecimal twoDecimal = BigDecimal.valueOf(2);
    protected BigInteger fore = BigInteger.valueOf(4);

    protected static final BigInteger RSA220 = new BigInteger("2260138526203405784941654048610197513508038915719776718321197768109445641817966676608593121306582577250631562886676970448070001811149711863002112487928199487482066070131066586646083327982803560379205391980139946496955261");
    protected static final BigInteger RSA1024 = new BigInteger("135066410865995223349603216278805969938881475605667027524485143851526510604859533833940287150571909441798207282164471551373680419703964191743046496589274256239341020864383202110372958725762358509643110564073501508187510676594629205563685529475213500852879416377328533906109750544334999811150056977236890927563");

    protected boolean logsEnabled = true;

    public void setLogsEnabled(boolean enabled) {
        this.logsEnabled = enabled;
    }

    private String emptyChar = (new char[1])[0] + "";

    protected void forceLog(Object... params){
        boolean logsEnabled = this.logsEnabled;
        this.logsEnabled = true;
        log(params);
        this.logsEnabled = logsEnabled;
    }

    protected void log(Object... params) {
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
        switch (millisecondsString.length()){
            case 1:
                millisecondsString = "00" + millisecondsString;
                break;
            case 2:
                millisecondsString = "0" + millisecondsString;
                break;
        }
        System.out.print(seconds + "." + millisecondsString + ": ");
    }

    protected void logInLine(Object... params) {
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
        int minLength = 10;
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
