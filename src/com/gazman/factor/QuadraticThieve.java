package com.gazman.factor;

import com.gazman.factor.matrix.BitMatrix;
import com.gazman.math.MathUtils;
import com.gazman.math.SqrRoot;

import java.math.BigInteger;
import java.util.ArrayList;

/**
 * Created by Ilya Gazman on 1/27/2016.
 */
public class QuadraticThieve extends Logger {
    private static final int B_SMOOTH = 4000;
    public static final double MINIMUM_LOG = 0.00000001;
    private int sieveVectorBound;
    private BigInteger primeBase[];
    private ArrayList<VectorData> vectorDatas = new ArrayList<>();
    private BigInteger N;
    private NumberData numberDatas[];
    private Wheel wheels[] = new Wheel[B_SMOOTH];
    private BigInteger root;

    public void factor(BigInteger input) {
        N = input;
        root = SqrRoot.bigIntSqRootCeil(input);

        log("Building Prime Base");
        buildPrimeBase();
        BigInteger highestPrime = primeBase[primeBase.length - 1];
//        baseLog = MathUtils.ln(new BigDecimal(root)).doubleValue();
        sieveVectorBound = highestPrime.intValue();
        log("Biggest prime is", highestPrime);
        log();

        log("Building sieve vector");
        initSieveVector();

        log("Building wheels");
        initSieveWheels();

        log("Start searching");
        long position = 0;
        long step = sieveVectorBound;

        while (true) {
//            double error = Math.log(primeBase[primeBase.length - 1].doubleValue());
            BigInteger target = root.add(BigInteger.valueOf(position)).pow(2).subtract(N);
            NumberData.upgrade();
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

    private void initSieveVector() {
        numberDatas = new NumberData[sieveVectorBound];
        for (int i = 0; i < numberDatas.length; i++) {
            numberDatas[i] = new NumberData(N, primeBase, root);
        }
    }

    private int sieve(long destination) {
        int newVectorsFound = 0;
        Double[] logs = new Double[sieveVectorBound];
//        BigInteger[] actualValues = new BigInteger[sieveVectorBound];
//        BigInteger[] originalValues = new BigInteger[sieveVectorBound];

        for (int i = 0; i < primeBase.length; i++) {
            Wheel wheel = wheels[i];
            wheel.savePosition();
            while (wheel.testMove(destination)) {
                double log = wheel.getLog();
//                MutableInteger mutableInteger = wheel.getMutableInteger();
                long position = wheel.move();
                int index = (int) (position % sieveVectorBound);
                if(logs[index] == null) {
                    BigInteger target = root.add(BigInteger.valueOf(position)).pow(2).subtract(N);
//                    actualValues[index] = target;
//                    originalValues[index] = target;
                    logs[index] = Math.log(target.doubleValue());
                }
                logs[index] -= log;

//                BigInteger divider = wheel.getDivider();
//                BigInteger actualValue = actualValues[index];
//                BigInteger originalValue = originalValues[index];
//                BigInteger target = root.add(BigInteger.valueOf(position)).pow(2).subtract(N);

//                log(originalValue, mutableInteger, target, "<", wheel, destination / sieveVectorBound, destination, index);

//                if(!new BigInteger(mutableInteger.toString()).equals(target)){
//                    throw new IllegalStateException("O.o " + wheel);
//                }
//                if(!actualValue.mod(divider).equals(BigInteger.ZERO)){
//                    try {
//                        Thread.sleep(100);
//                    } catch (InterruptedException e) {
//                        e.printStackTrace();
//                    }
//                    throw new IllegalStateException("O.o " + wheel + " " + divider + " " + actualValue);
//                }
//                actualValues[index] = actualValues[index].divide(wheel.getDivider());
            }
            wheel.restorePosition();
        }

//        int expected = 0;
//        for (int i = 0; i < logs.length; i++) {
//            if (logs[i] != null && logs[i] < 0.00000001) {
//                expected++;
//            }
//        }
//        log("Expected", expected);
//        expected = expected;

        for (int i = 0; i < primeBase.length; i++) {
            Wheel wheel = wheels[i];

            while (wheel.testMove(destination)) {
                long position = wheel.move();
                int index = (int) (position % sieveVectorBound);

                if (logs[index] > MINIMUM_LOG) {
                    continue;
                }

                int powers = numberDatas[index].sieve(i, position);
                if (powers > 0) {
                    vectorDatas.add(numberDatas[index].getData());
                    newVectorsFound++;
                }
            }
        }

//        if (newVectorsFound > 0) {
//            log("Expected", expected, "Found", newVectorsFound, "total", vectorDatas.size());
//        }


        return newVectorsFound;
    }

    private boolean tryToSolve() {
        if (vectorDatas.size() < B_SMOOTH) {
            log("Found", vectorDatas.size());
            return false;
        }

        log("Building vectors");
        for (VectorData vectorData : vectorDatas) {
            vectorData.buildVector(primeBase);
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
            x = x.multiply(vectorData.x).mod(N);
            y = y.multiply(vectorData.y);
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
