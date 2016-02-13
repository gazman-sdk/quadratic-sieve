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
    private static final int B_SMOOTH = 10000;
    private static final double MINIMUM_LOG = 0.0000001;
    private  double MINIMUM_ESTIMATED_LOG = 2;
    private int sieveVectorBound;
    private BigInteger primeBase[];
    private ArrayList<VectorData> vectorDatas = new ArrayList<>();
    private BigInteger N;
    private Wheel wheels[] = new Wheel[B_SMOOTH];
    private BigInteger root;
    private double baseLog;

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
            baseLog = calculateBaseLog(position);
            position += step;
            MINIMUM_ESTIMATED_LOG = calculateBaseLog(position) - baseLog;
            if (sieve(position) > 0 && tryToSolve()) {
                break;
            }
        }
    }

    private double calculateBaseLog(long position) {
        double target = root.add(BigInteger.valueOf(position)).pow(2).subtract(N).doubleValue();
        return Math.log(target);
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
        Double[] trueLogs = new Double[sieveVectorBound];
        BitSet[] vectors = new BitSet[sieveVectorBound];

        for (int i = 0; i < primeBase.length; i++) {
            Wheel wheel = wheels[i];
            wheel.savePosition();
            while (wheel.testMove(destination)) {
                double log = wheel.nextLog();
                long position = wheel.move();
                int index = (int) (position % sieveVectorBound);
                if (logs[index] == null) {
                    logs[index] = 0d;
                }
                logs[index] += log;
            }
            wheel.restorePosition();
        }

        for (int i = primeBase.length - 1; i >= 0; i--) {
            Wheel wheel = wheels[i];
            while (wheel.testMove(destination)) {
                wheel.nextLog();
                long position = wheel.move();
                int index = (int) (position % sieveVectorBound);
                if (logs[index] == null){
                    continue;
                }
                if (trueLogs[index] == null) {
                    if (baseLog - logs[index] > MINIMUM_ESTIMATED_LOG) {
                        continue;
                    }
                    double trueLog = calculateBaseLog(position);
                    trueLogs[index] = trueLog;
                }

                if (trueLogs[index] - logs[index] > MINIMUM_LOG) {
                    continue;
                }

                if (vectors[index] == null) {
                    vectors[index] = new BitSet(i);
                    vectorDatas.add(new VectorData(vectors[index], index + destination - sieveVectorBound));
                    newVectorsFound++;
                }
                if (wheel.getPowers() % 2 != 0) {
                    vectors[index].set(i);
                }
            }
        }

        return newVectorsFound;
    }

    private boolean tryToSolve() {
        log("Found", vectorDatas.size());
        if (vectorDatas.size() < B_SMOOTH) {
            return false;
        }

        BitMatrix bitMatrix = new BitMatrix();
        ArrayList<ArrayList<VectorData>> solutions = bitMatrix.solve(vectorDatas);

        for (int i = 0; i < solutions.size(); i++) {
            ArrayList<VectorData> solution = solutions.get(i);
            log("Testing solution", (i + 1) + "/" + solutions.size());
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

        for (int i = 0; i < B_SMOOTH; ) {
            prime = prime.nextProbablePrime();
            if (MathUtils.isRootInQuadraticResidues(N, prime)) {
                primeBase[i] = prime;
                i++;
            }
        }
    }
}
