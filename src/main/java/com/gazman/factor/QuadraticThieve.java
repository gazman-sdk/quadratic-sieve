package com.gazman.factor;

import com.gazman.factor.matrix.BitMatrix;
import com.gazman.factor.matrix.VectorsShrinker;
import com.gazman.factor.wheels.Wheel;
import com.gazman.math.MathUtils;
import com.gazman.math.SqrRoot;

import java.math.BigInteger;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Created by Ilya Gazman on 1/27/2016.
 */
public class QuadraticThieve extends Logger {
    private static final int B_SMOOTH = 5000;
    private static final double MINIMUM_LOG = 0.0000001;
    public static final int MAX_LOOPS = B_SMOOTH * 2;
    public static final int LOGS_TIME_BY_LOOPS = B_SMOOTH / 20;

    private final double minimumBigPrimeLog;
    private final int sieveVectorBound;
    private final BigInteger primeBase[] = new BigInteger[B_SMOOTH];
    private final int step;
    private ArrayList<VectorData> bSmoothVectors = new ArrayList<>();
    private final BigInteger N;
    private final BigInteger root;
    private int bSmoothFound;
    private BigPrimesList bigPrimesList = new BigPrimesList();
    private VectorsShrinker vectorsShrinker = new VectorsShrinker();
    private final double double2Root;
//    private final int threadCount = Runtime.getRuntime().availableProcessors(); not working on mac
    private final int threadCount = 2;
    private final AtomicInteger speedCounter = new AtomicInteger(0);
    private final AtomicInteger speed = new AtomicInteger(0);
    private long startingTime;

    public QuadraticThieve(BigInteger input) {
        log("Factoring started");
        N = input;
        root = SqrRoot.bigIntSqRootCeil(input);
        double2Root = root.add(root).doubleValue();

        log("Building Prime Base");
        buildPrimeBase();
        vectorsShrinker.init(root, primeBase.length, N);
        BigInteger highestPrime = primeBase[primeBase.length - 1];
        sieveVectorBound = highestPrime.intValue();
        minimumBigPrimeLog = Math.log(highestPrime.pow(2).doubleValue());
        step = sieveVectorBound;
        log("Biggest prime is", highestPrime);
        log();

        log("Working on", threadCount, "threads");
        log("Start searching");
    }

    public void start() {
        for (int i = 0; i < threadCount; i++) {
            int threadId = i;
            new Thread(() -> execute(threadId)).start();
        }
    }

    private void execute(int threadId) {
        long basePosition = 0;
        if(threadId == 0){
            startingTime = System.currentTimeMillis();
        }
        while (true) {
            long position = basePosition + threadId * MAX_LOOPS * step;
            log(threadId, "Building wheels");
            Wheel[] localWheels = initSieveWheels(position);
            log(threadId, "Started");
            for (int loops = 0; loops < MAX_LOOPS; loops++) {
                double baseLog = calculateBaseLog(position);
                position += step;
                boolean sieve = sieve(position, baseLog, localWheels);
                speed.incrementAndGet();
                if(speedCounter.incrementAndGet() == LOGS_TIME_BY_LOOPS){
                    speedCounter.set(0);
                    logProcesses();
                }

                if (sieve && isReadyToBeSolved()) {
                    log(threadId, "Getting ready to solve");
                    if (threadId == 0) {
                        try {
                            Thread.sleep(10);
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        }
                    } else {
                        synchronized (QuadraticThieve.this) {
                            try {
                                log(threadId, "Waiting for solution");
                                QuadraticThieve.this.wait();
                                continue;
                            } catch (InterruptedException e) {
                                e.printStackTrace();
                            }
                        }
                    }
                    logProcesses();
                    if (tryToSolve()) {
                        System.exit(1);
                        break;
                    } else {
                        System.exit(1);
                        synchronized (QuadraticThieve.this) {
                            QuadraticThieve.this.notifyAll();
                        }
                    }
                }
            }
            basePosition += step * MAX_LOOPS * threadCount;
        }
    }

    private void logProcesses() {
        int speed = this.speed.intValue();
        long currentTimeMillis = System.currentTimeMillis();
        long secPass = (currentTimeMillis - startingTime) / 1000;
        if(secPass == 0){
            return;
        }
        log("speed", speed / secPass * step / 1000, "kValues a second, B-Smooth found", bSmoothFound, "Big primes found", bigPrimesList.getPrimesFound());
    }

    private double calculateBaseLog(double position) {
//        double target = root.add(BigInteger.valueOf(position)).pow(2).subtract(N).doubleValue();
        return Math.log(position * (position + double2Root));

    }

    private Wheel[] initSieveWheels(long position) {
        Wheel[] wheels = new Wheel[B_SMOOTH];
        for (int i = 0; i < wheels.length; i++) {
            wheels[i] = new Wheel(primeBase[i], N, root.add(BigInteger.valueOf(position)), sieveVectorBound);
        }
        return wheels;
    }

    private boolean sieve(long destination, double baseLog, Wheel[] wheels) {
        boolean vectorsFound = false;
        double[] logs = new double[sieveVectorBound];
        double[] trueLogs = new double[sieveVectorBound];
        VectorData[] vectors = new VectorData[sieveVectorBound];

        for (int i = 0; i < primeBase.length; i++) {
            Wheel wheel = wheels[i];
            wheel.savePosition();
            wheel.prepareToMove();
            while (wheel.testMove()) {
                int index = wheel.move();
                if(index > logs.length){
                    log(index - sieveVectorBound, logs.length - sieveVectorBound, "error");
                    System.exit(3);
                }
                logs[index] += wheel.log;
            }
            wheel.restorePosition();
        }

        for (int i = primeBase.length - 1; i >= 0; i--) {
            Wheel wheel = wheels[i];
            wheel.prepareToMove();
            while (wheel.testMove()) {
                int index = wheel.move();

                if (trueLogs[index] == 0) {
                    if (baseLog - logs[index] > minimumBigPrimeLog) {
                        continue;
                    }
                    trueLogs[index] = calculateBaseLog(destination + index - sieveVectorBound);
                }

                double reminderLog = trueLogs[index] - logs[index];
                if (reminderLog > minimumBigPrimeLog) {
                    continue;
                }

                boolean bigPrime = reminderLog > MINIMUM_LOG;

                synchronized (this) {
                    if (vectors[index] == null) {
                        VectorData vectorData = new VectorData(new BitSet(i), index + destination - sieveVectorBound);
                        vectors[index] = vectorData;

                        if (bigPrime) {
                            long prime = Math.round(Math.pow(Math.E, reminderLog));
                            bigPrimesList.add(prime, vectorData);
                        } else {
                            bSmoothVectors.add(vectorData);
                            bSmoothFound++;
                        }
                    }
                    vectorsFound = true;
                    vectors[index].vector.set(i);
                }
            }
        }

        return vectorsFound;
    }

    private boolean tryToSolve() {
        log("Building matrix");

        ArrayList<VectorData> vectorDatas = vectorsShrinker.shrink(bSmoothVectors, bigPrimesList);

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

    private boolean isReadyToBeSolved() {
        return bSmoothVectors.size() + bigPrimesList.getPrimesFound() >= B_SMOOTH;
    }

    private boolean testSolution(ArrayList<VectorData> solutionVector) {
        BigInteger y = one;
        BigInteger x = one;

        for (VectorData vectorData : solutionVector) {
            BigInteger savedX, savedY;
            if (vectorData.x != null) {
                savedX = vectorData.x;
                savedY = vectorData.y;
            } else {
                savedX = root.add(BigInteger.valueOf(vectorData.position));
                savedY = savedX.pow(2).subtract(N);
            }
            x = x.multiply(savedX).mod(N);
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
