package com.gazman.factor;

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

        int solutions = 0;
        int solutionA = -1, solutionB = -1;
        for (int i = 0; i < this.prime; i++) {
            if ((i * i) % this.prime == target) {
                solutions++;
                switch (solutions) {
                    case 1:
                        solutionA = i;
                        break;
                    case 2:
                        solutionB = i;
                        break;
                }
                if (solutions == 2) {
                    break;
                }
            }
        }

        positions = new long[solutions];
        if(solutions > 0){
            positions[0] = solutionA - root.mod(prime).intValue();
        }
        if (solutions > 1) {
            positions[1] = solutionB - root.mod(prime).intValue();
        }
        for (int i = 0; i < solutions; i++) {
            if(positions[i] < 0){
                positions[i] += this.prime;
            }
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
