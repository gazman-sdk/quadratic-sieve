package com.gazman.factor.wheels;

import com.gazman.math.MathUtils;

import java.math.BigInteger;

/**
 * Created by Ilya Gazman on 2/3/2016.
 */
public class Wheel implements Cloneable {
    private int[] positions;
    private int count = 0;
    private int prime;
    public final double log;
    private int[] savedPosition;
    private int savedCount;
    private int sieveVectorBound;
    private int[] loops;

    public Wheel(BigInteger prime, BigInteger N, BigInteger root, int sieveVectorBound) {
        this.sieveVectorBound = sieveVectorBound;
        int target = N.mod(prime).intValue();
        this.prime = prime.intValue();
        log = Math.log(this.prime);

        long[] flats = MathUtils.ressol(this.prime, target);

        int solutions = 0;
        int flatsLength = flats.length;
        //noinspection ForLoopReplaceableByForEach
        for (int i = 0; i < flatsLength; i++) {
            long flat = flats[i];
            if (flat > -1) {
                solutions++;
            }
        }

        positions = new int[solutions];
        savedPosition = new int[solutions];
        loops = new int[solutions];
        for (int i = 0; i < flatsLength; i++) {
            if (flats[i] <= -1) {
                continue;
            }
            positions[i] = (int) (flats[i] - root.mod(prime).intValue());
            if (positions[i] < 0) {
                positions[i] += this.prime;
            }
        }

        if (positions.length == 2 && positions[0] > positions[1]) {
            int tmp = positions[0];
            positions[0] = positions[1];
            positions[1] = tmp;
        }
    }

    public int move() {
        int position = positions[count];
        positions[count] += prime;
        if (positions[count] >= sieveVectorBound) {
            positions[count] -= sieveVectorBound;
            loops[count]--;
        }
        count++;
        if (count >= positions.length) {
            count -= positions.length;
        }

        return position;
    }

    public void seek(long destination) {
        if (positions.length == 1) {
            positions[0] += (prime * destination);
        } else {
            positions[0] += destination / 2 * prime;
            positions[1] += destination / 2 * prime;
        }
        if (destination % 2 == 1) {
            positions[0] += prime;
            count++;
            if (count >= positions.length) {
                count -= positions.length;
            }
        }
        for (int i = 0; i < positions.length; i++) {
            positions[i] = positions[i] % sieveVectorBound;
        }
    }

    public boolean testMove() {
        return loops[count] > 0;
    }

    public void savePosition() {
        savedPosition[0] = positions[0];
        if (positions.length > 1) {
            savedPosition[1] = positions[1];
        }
        savedCount = count;
    }

    public void restorePosition() {
        positions[0] = savedPosition[0];
        if (positions.length > 1) {
            positions[1] = savedPosition[1];
        }
        count = savedCount;
    }

    @Override
    public String toString() {
        return prime + "";
    }

    public void prepareToMove() {
        for (int i = 0; i < loops.length; i++) {
            loops[i]++;
        }
    }

    @Override
    public Wheel clone() {
        try {
            return (Wheel) super.clone();
        } catch (CloneNotSupportedException e) {
            e.printStackTrace();
        }
        return null;
    }
}
