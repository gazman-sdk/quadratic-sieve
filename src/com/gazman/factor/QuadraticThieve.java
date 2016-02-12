package com.gazman.factor;

import com.gazman.factor.matrix.BitMatrix;
import com.gazman.math.MathUtils;
import com.gazman.math.SqrRoot;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.BitSet;

/**
 * Created by Ilya Gazman on 1/27/2016.
 */
public class QuadraticThieve extends Logger {
    private static final int B_SMOOTH = 100000;
    public static final double MINIMUM_LOG = 0.00000001;
    private int sieveVectorBound;
    private BigInteger primeBase[];
    private ArrayList<VectorData> vectorDatas = new ArrayList<>();
    private BigInteger N;
    private Wheel wheels[] = new Wheel[B_SMOOTH];
    private BigInteger root;

    public void factor(BigInteger input) {
        log("Factoring started");
        N = input;
        root = SqrRoot.bigIntSqRootCeil(input);

        log("Building Prime Base");
        buildPrimeBase();
        BigInteger highestPrime = primeBase[primeBase.length - 1];
        sieveVectorBound = highestPrime.intValue();
        log("Biggest prime is", highestPrime);
        log();

        log("Building wheels");
        initSieveWheels();

        log("Start searching");
        long position = 0;
        long step = sieveVectorBound;

        while (true) {
            position += step;
            if (sieve(position) > 0 && tryToSolve()) {
                break;
            }
        }
    }

    private void initSieveWheels() {
        for (int i = 0; i < wheels.length; i++) {
            wheels[i] = new Wheel();
            wheels[i].init(primeBase[i], N, root);
        }
    }

    private int sieve(long destination) {
        int newVectorsFound = 0;
        Double[] logs = new Double[sieveVectorBound];
        BitSet[] vectors = new BitSet[sieveVectorBound];

        for (int i = 0; i < primeBase.length; i++) {
            Wheel wheel = wheels[i];
//            wheel.savePosition();
            while (wheel.testMove(destination)) {
                double log = wheel.getLog();
                long position = wheel.move();
                int index = (int) (position % sieveVectorBound);
                if(logs[index] == null) {
                    BigInteger target = root.add(BigInteger.valueOf(position)).pow(2).subtract(N);
                    logs[index] = Math.log(target.doubleValue());
                    vectors[index] = new BitSet();
                }
                logs[index] -= log;
                if(wheel.getPowers() % 2 != 0){
                    vectors[index].set(i);
                }
            }
//            wheel.restorePosition();
        }

        for (int i = 0; i < primeBase.length; i++) {
            if (logs[i] != null && logs[i] < MINIMUM_LOG) {
                vectorDatas.add(new VectorData(vectors[i], i + destination - sieveVectorBound));
                newVectorsFound++;
            }
        }

        return newVectorsFound;
    }

    private boolean tryToSolve() {
        if (vectorDatas.size() < B_SMOOTH) {
            log("Found", vectorDatas.size());
            return false;
        }

        log("Solving...", vectorDatas.size());

        BitMatrix bitMatrix = new BitMatrix();
//        bitMatrix.setLogsEnabled(false);
        ArrayList<ArrayList<VectorData>> solutions2 = bitMatrix.solve(vectorDatas);

        for (ArrayList<VectorData> solution : solutions2) {
            log("Testing solution");
            if (testSolution(solution)) {
                return true;
            }
        }
        log("no luck");

        return false;
    }

    private boolean testSolution(ArrayList<VectorData> solutionVector) {
        BigInteger y = one;
        BigInteger x = one;

        for (VectorData vectorData : solutionVector) {
            BigInteger savedX = root.add(BigInteger.valueOf(vectorData.position));
            BigInteger savedY = savedX.pow(2).subtract(N);
            x = x.multiply(savedX);
            y = y.multiply(savedY);
        }

        y = SqrRoot.bigIntSqRootFloor(y);
        BigInteger gcd = N.gcd(x.add(y));
        if (!gcd.equals(one) && !gcd.equals(N)) {
            log("Solved");
            log(gcd);

            return true;
        }

        return false;
    }

    private void buildPrimeBase() {
        primeBase = new BigInteger[B_SMOOTH];
        BigInteger prime = BigInteger.ONE;

//        prime = prime.nextProbablePrime();

        for (int i = 0; i < B_SMOOTH; ) {
            prime = prime.nextProbablePrime();
            if (MathUtils.isRootInQuadraticResidues(N, prime)) {
                primeBase[i] = prime;
                i++;
            }
        }
    }
}
