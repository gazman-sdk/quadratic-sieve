package com.gazman.factor.wheels;

import com.gazman.factor.MutableInteger;

import java.math.BigInteger;

/**
 * Created by Ilya Gazman on 2/3/2016.
 */
public class PowerWheel extends Wheel {
    private MutableInteger mutableInteger[] = new MutableInteger[2];
    private MutableInteger savedMutableInteger[];
    private int powers = 1;

    public void init(BigInteger prime, BigInteger N, BigInteger root) {
        super.init(prime, N, root);

        for (int i = 0; i < positions.length; i++) {
            BigInteger bigPosition = BigInteger.valueOf(positions[i]);
            mutableInteger[i] = new MutableInteger(this.prime, root.add(bigPosition), N);
        }
    }

    public double nextLog() {
        powers = mutableInteger[count].nextPower();
        return super.nextLog() * powers;
    }

    public int getPowers() {
        return powers;
    }

    public void savePosition() {
        super.savePosition();
        savedMutableInteger = new MutableInteger[mutableInteger.length];
        for (int i = 0; i < mutableInteger.length; i++) {
            if (mutableInteger[i] != null) {
                savedMutableInteger[i] = mutableInteger[i].clone();
            }
        }
    }

    public void restorePosition() {
        super.restorePosition();
        mutableInteger = savedMutableInteger;
    }
}
