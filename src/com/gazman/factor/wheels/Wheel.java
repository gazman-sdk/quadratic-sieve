package com.gazman.factor.wheels;

import com.gazman.math.MathUtils;

import java.math.BigInteger;

/**
 * Created by Ilya Gazman on 2/3/2016.
 */
public class Wheel {
    protected long[] positions;
    protected int count = 0;
    protected int prime;
    private double log;
    private long[] savedPosition;
    private int savedCount;

    public void init(BigInteger prime, BigInteger N, BigInteger root) {
        int target = N.mod(prime).intValue();
        this.prime = prime.intValue();
        log = Math.log(this.prime);

        long[] flats = MathUtils.ressol(this.prime, target);

        int solutions =  0;
        for (int i = 0; i < flats.length; i++) {
            if (flats[i] > -1) {
                solutions++;
            }
        }

        positions = new long[solutions];
        for (int i = 0; i < flats.length; i++) {
            if (flats[i] <= -1) {
                continue;
            }
            positions[i] = flats[i] - root.mod(prime).intValue();
            if (positions[i] < 0) {
                positions[i] += this.prime;
            }
        }

        if (positions.length == 2 && positions[0] > positions[1]) {
            long tmp = positions[0];
            positions[0] = positions[1];
            positions[1] = tmp;
        }
    }

    public double nextLog() {
        return log;
    }

    public int getPowers() {
        return 1;
    }

    public long move() {
        long position = positions[count];
        positions[count] += prime;
        count = (count + 1) % positions.length;
        return position;
    }

    public boolean testMove(long limit) {
        return positions.length != 0 && positions[count] < limit;
    }

    public void savePosition() {
        savedPosition = positions.clone();
        savedCount = count;
    }

    public void restorePosition() {
        positions = savedPosition.clone();
        count = savedCount;
    }

    @Override
    public String toString() {
        return prime + "";
    }
}
