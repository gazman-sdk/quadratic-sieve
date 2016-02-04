package com.gazman.factor;

import com.gazman.math.MathUtils;
import com.gazman.math.SqrRoot;

import java.math.BigInteger;
import java.util.ArrayList;

/**
 * Created by Ilya Gazman on 1/27/2016.
 */
public class QuadraticThieve extends BaseFactor {
    private static final int B_SMOOTH = 50000;
    private static final int SIEVE_VECTOR_BOUND = 1000000;
    private BigInteger primeBase[];
    private ArrayList<VectorData> vectorDatas = new ArrayList<>();
    private BigInteger N;
    private Matrix matrix = new Matrix();
    private VectorNumber vectorNumbers[] = new VectorNumber[SIEVE_VECTOR_BOUND];
    private Wheel wheels[] = new Wheel[B_SMOOTH];
    private BigInteger root;

    public void factor(BigInteger input) {
        N = input;
        root = SqrRoot.bigIntSqRootCeil(input);
        long position = 0;
        long step =  SIEVE_VECTOR_BOUND;

        matrix.setLogsEnabled(false);

        log("Building Prime Base");
        buildPrimeBase();
        log("Biggest prime is", primeBase[primeBase.length - 1]);
        log("Building sieve vector");
        initSieveVector();
        log("Building wheels");
        initSieveWheels();

        log("Start searching");

        while (true) {
            VectorNumber.upgrade();
            position += step;
            if (sieve(position) && tryToSolve()) {
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
        for (int i = 0; i < vectorNumbers.length; i++) {
            vectorNumbers[i] = new VectorNumber(N, primeBase, root);
        }
    }

    private boolean sieve(long destination) {
        boolean foundNewVectors = false;
        for (int i = 0; i < primeBase.length; i++) {
            Wheel wheel = wheels[i];

            while (wheel.testMove(destination)){
                long position = wheel.move();
                int index = (int) (position % SIEVE_VECTOR_BOUND);
                if(vectorNumbers[index].sieve(i, position)){
                    vectorDatas.add(vectorNumbers[index].getData());
                    foundNewVectors = true;
                }
            }
        }
        return foundNewVectors;
    }

    private boolean tryToSolve() {
        if (vectorDatas.size() < B_SMOOTH) {
            log("Found vector", vectorDatas.size());
            return false;
        }

        for (VectorData vectorData : vectorDatas) {
            vectorData.buildVector(primeBase);
        }

        log("Solving...", vectorDatas.size());
        matrix = new Matrix();
        matrix.setLogsEnabled(false);
        ArrayList<Object> solutions = matrix.solve(vectorDatas);
        for (Object solution : solutions) {
            if (testSolution((boolean[]) solution)) {
                return true;
            }
        }
        log("no luck");

        return false;
    }

    private boolean testSolution(boolean[] solutionVector) {
        BigInteger y = one;
        BigInteger x = one;
        for (int i = 0; i < solutionVector.length; i++) {
            if (solutionVector[i]) {
                VectorData vectorData = vectorDatas.get(i);
                y = y.multiply(vectorData.y);
                x = x.multiply(vectorData.x);
            }
        }

        //noinspection SuspiciousNameCombination
        y = SqrRoot.bigIntSqRootFloor(y);
        x = SqrRoot.bigIntSqRootFloor(x);
        BigInteger gcd = N.gcd(x.subtract(y));
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
