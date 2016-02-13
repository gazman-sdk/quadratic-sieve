package com.gazman.factor;

import com.gazman.math.MathUtils;

import java.math.BigInteger;

/**
 * Created by Ilya Gazman on 2/9/2016.
 */
public class MutableInteger implements Cloneable{
    private int[] digits;
    private int[] delta;
    private final int[] deltaDelta;
    private static final BigInteger TOW = BigInteger.valueOf(2);
    private int base;

    private MutableInteger(int[] deltaDelta){
        this.deltaDelta = deltaDelta;
    }

    public MutableInteger(int base, BigInteger root, BigInteger N) {
        this.base = base;
        this.digits = toLocalBase(root.pow(2).subtract(N));
        BigInteger bigBase = BigInteger.valueOf(base);
        delta = toLocalBase(root.multiply(TOW).multiply(bigBase).add(bigBase.pow(2)));
        deltaDelta = toLocalBase(bigBase.pow(2).multiply(BigInteger.valueOf(2)));
        if(digits.length < delta.length){
            moreDigits.update();
        }
    }

    private SizeWatcher moreDigits = new SizeWatcher() {

        @Override
        protected int[] update() {
            digits = doubleSize(digits);
            return digits;
        }
    };
    private SizeWatcher moreDelta = new SizeWatcher() {
        @Override
        protected int[] update() {
            delta = doubleSize(delta);
            return delta;
        }
    };

    private abstract class SizeWatcher {

        protected abstract int[] update();

        protected int[] doubleSize(int[] digits) {
            int[] newDigits = new int[digits.length * 2];
            System.arraycopy(digits, 0, newDigits, 0, digits.length);
            return newDigits;
        }
    }

    public int nextPower() {
        add(digits, delta, moreDigits);
        add(delta, deltaDelta, moreDelta);
        int i = 0;
        //noinspection StatementWithEmptyBody
        for (; i < digits.length && digits[i] == 0; i++)
            ;
        return i;
    }

    private void add(int[] digits, int[] delta, SizeWatcher sizeWatcher) {
        int extra = 0;
        for (int i = 0; (i < delta.length || extra == 1); i++) {
            if(i == digits.length){
                digits = sizeWatcher.update();
            }
            digits[i] = digits[i] + (i < delta.length ? delta[i] : 0) + extra;
            if (digits[i] >= base) {
                digits[i] = digits[i] - base;
                extra = 1;
            } else {
                extra = 0;
            }
        }
    }

    private int[] toLocalBase(BigInteger value) {
        int[] digits = new int[(int) (Math.round(MathUtils.log(value.doubleValue(), base)) + 1)];
        BigInteger bigCycle = BigInteger.valueOf(base);
        for (int i = 0; i < digits.length; i++) {
            BigInteger mod = value.mod(bigCycle);
            int intMod = mod.intValue();
            digits[i] = intMod;
            value = value.subtract(mod).divide(bigCycle);
        }
        return digits;
    }

    @Override
    public String toString() {
        return digitsToString(this.digits);
    }

    private String digitsToString(int[] digits) {

        BigInteger result = BigInteger.ZERO;
        BigInteger base = BigInteger.valueOf(this.base);
        for (int i = 0; i < digits.length; i++) {
            result = result.add(base.pow(i).multiply(BigInteger.valueOf(digits[i])));
        }

        return result.toString();
    }

    @SuppressWarnings("CloneDoesntCallSuperClone")
    @Override
    public MutableInteger clone() {
        MutableInteger mutableInteger = new MutableInteger(deltaDelta);
        mutableInteger.digits = digits.clone();
        mutableInteger.delta = delta.clone();
        mutableInteger.base = base;
        return mutableInteger;
    }
}
