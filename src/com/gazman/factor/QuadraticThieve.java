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
    private static final int B_SMOOTH = 1000;
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
        sieveVectorBound = primeBase[primeBase.length - 1].intValue();
        log("Biggest prime is", primeBase[primeBase.length - 1]);
        log();

        log("Building sieve vector");
        initSieveVector();

        log("Building wheels");
        initSieveWheels();

        log("Start searching");
        long position = 0;
        long step = sieveVectorBound;

        while (true) {
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
        int foundNewVectors = 0;
        int maxDivisionsIndex = -1;
        int maxDivisions = 0;
        for (int i = 0; i < primeBase.length; i++) {
            Wheel wheel = wheels[i];

            while (wheel.testMove(destination)) {
                long position = wheel.move();
                int index = (int) (position % sieveVectorBound);
                if (numberDatas[index].sieve(i, position)) {
                    vectorDatas.add(numberDatas[index].getData());
                    foundNewVectors++;
                } else if (numberDatas[index].getDivisions() > maxDivisions) {
                    maxDivisions = numberDatas[index].getDivisions();
                    maxDivisionsIndex = index;
                }
            }
        }
        return foundNewVectors > 0 ? foundNewVectors : -maxDivisionsIndex;
    }

    private boolean tryToSolve() {
        if (vectorDatas.size() < B_SMOOTH) {
            log("Found vector", vectorDatas.size());
            return false;
        }

        log("Building vectors");
        for (VectorData vectorData : vectorDatas) {
            vectorData.buildVector(primeBase);
        }

        log("Solving...", vectorDatas.size());

        BitMatrix bitMatrix = new BitMatrix();
        bitMatrix.setLogsEnabled(false);
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
            x = x.multiply(vectorData.x);
            y = y.multiply(vectorData.y);
        }

        y = SqrRoot.bigIntSqRootFloor(y);
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
