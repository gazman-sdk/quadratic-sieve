package com.gazman.factor;

import com.gazman.math.MathUtils;
import com.gazman.math.SqrRoot;

import java.math.BigInteger;
import java.util.ArrayList;

/**
 * Created by Ilya Gazman on 1/27/2016.
 */
public class QuadraticThieve extends BaseFactor {
    private static final int B_SMOOTH = 1000;
    private int sieveVectorBound;
    private BigInteger primeBase[];
    private ArrayList<VectorData> vectorDatas = new ArrayList<>();
    private BigInteger N;
    private Matrix matrix = new Matrix();
    private NumberData numberDatas[];
    private Wheel wheels[] = new Wheel[B_SMOOTH];
    private BigInteger root;

    public void factor(BigInteger input) {
        N = input;
        root = SqrRoot.bigIntSqRootCeil(input);
        matrix.setLogsEnabled(false);

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
            int newVectorsFromSieving = sieve(position);
            if (newVectorsFromSieving < 0) {
//                searchForGiants(position, newVectorsFromSieving * -1);
            } else {
                if (tryToSolve()) {
                    break;
                }
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

    private void searchForGiants(long position, int giantIndex) {
        BigInteger giant = numberDatas[giantIndex].getRemindY();
        BigInteger nextGiant = zero;
        BigInteger x = root.
                add(BigInteger.valueOf(position + giantIndex - sieveVectorBound)).
                add(giant).
                pow(2);
        nextGiant = x.subtract(N);
        if (!nextGiant.mod(giant).equals(zero)) {
            throw new IllegalStateException("O.o");
        }

        nextGiant = nextGiant.divide(giant);

        for (int i = 0; i < primeBase.length; i++) {
            BigInteger prime = primeBase[i];
            while (nextGiant.mod(prime).equals(zero)){
                nextGiant = nextGiant.divide(prime);
            }
            if(nextGiant.equals(one)){
                log("Found giant");

                break;
            }
        }
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
