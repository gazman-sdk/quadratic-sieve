package com.gazman.factor;

import com.gazman.math.MathUtils;

import java.math.BigInteger;

/**
 * Created by Ilya Gazman on 2/3/2016.
 */
public class Wheel {
    private long[] positions;
    private int count = 0;
    private int prime;

    public void init(BigInteger prime, BigInteger N, BigInteger root) {
        int target = N.mod(prime).intValue();
        this.prime = prime.intValue();


        long[] flats = MathUtils.ressol(this.prime, target);

        int solutions = 0;
        for (int i = 0; i < flats.length; i++) {
            if(flats[i] > -1){
                solutions++;
            }
        }

        positions = new long[solutions];
        for (int i = 0; i < flats.length; i++) {
            if(flats[i] > -1){
                positions[i] = flats[i] - root.mod(prime).intValue();
            }
        }
        for (int i = 0; i < solutions; i++) {
            if(positions[i] < 0){
                positions[i] += this.prime;
            }
        }
        if(positions.length == 2 && positions[0] > positions[1]){
            long tmp = positions[0];
            positions[0] = positions[1];
            positions[1] = tmp;
        }
    }

    public long move() {
        long position = positions[count];
        positions[count] += prime;
        count = (count + 1) % positions.length;
        return position;
    }

    public boolean testMove(long limit) {
        if(positions.length ==0){
            return false;
        }
        return positions[count] < limit;
    }

}
