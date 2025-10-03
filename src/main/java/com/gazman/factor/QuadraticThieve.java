package com.gazman.factor;

import com.gazman.factor.matrix.BitMatrix;
import com.gazman.factor.matrix.VectorsShrinker;
import com.gazman.factor.wheels.Wheel;
import com.gazman.math.MathUtils;
import com.gazman.math.SqrRoot;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Quadratic Sieve with:
 *  - Single-pass sieve (scaled logs)
 *  - Per-index bucket lists of hitting primes (massive TD reduction)
 *  - Single-large-prime policy aligned to original (<= B_max^2)
 *  - Per-thread batching (no hot-path sync/alloc)
 */
public class QuadraticThieve extends Logger {
    // --- Tuning knobs ---
    private static final int B_SMOOTH = 5000;               // factor base size
    public static final int MAX_LOOPS = B_SMOOTH * 2;       // blocks per thread before moving window
    public static final int LOGS_TIME_BY_LOOPS = B_SMOOTH / 20;

    // Fixed-point logs for hot loop
    private static final int LOG_SCALE = 256;
    private static final int LOG_EPS = 6;                   // ~0.023 nats

    // Growth factor for bucket storage
    private static final int BUCKET_GROWTH = 2;

    private final BigInteger N;
    private final BigInteger root;
    private final double double2Root;

    private final BigInteger[] primeBase = new BigInteger[B_SMOOTH];
    private final int[] primeBaseInt = new int[B_SMOOTH];
    private final int[] logPScaled = new int[B_SMOOTH];

    private final int sieveVectorBound; // block size
    private final int step;

    // Large-prime cutoff = (B_max)^2 (like original)
    private final BigInteger largePrimeBoundBI;
    private final int BIG_PRIME_CUTOFF_SCALED;

    // Threads
    private final int threadCount = 2;
    private final VectorsShrinker vectorsShrinker = new VectorsShrinker();
    private final BigPrimesList bigPrimesList = new BigPrimesList();

    private final ArrayList<VectorData> bSmoothVectors = new ArrayList<>();

    private final AtomicInteger speedCounter = new AtomicInteger(0);
    private final AtomicInteger speed = new AtomicInteger(0);
    private volatile long startingTime;

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
        step = sieveVectorBound;

        for (int i = 0; i < B_SMOOTH; i++) {
            int p = primeBase[i].intValue();
            primeBaseInt[i] = p;
            logPScaled[i] = (int) Math.round(Math.log(p) * LOG_SCALE);
        }

        largePrimeBoundBI = highestPrime.multiply(highestPrime);
        BIG_PRIME_CUTOFF_SCALED = (int) Math.round(Math.log(largePrimeBoundBI.doubleValue()) * LOG_SCALE);

        log("Biggest prime is", highestPrime);
        log();
        log("Working on", threadCount, "threads");
        log("Start searching");
    }

    public void start() {
        for (int i = 0; i < threadCount; i++) {
            final int threadId = i;
            new Thread(() -> execute(threadId), "QS-" + threadId).start();
        }
    }

    private void execute(int threadId) {
        long basePosition = 0;
        if (threadId == 0) startingTime = System.currentTimeMillis();

        while (true) {
            long position = basePosition + (long) threadId * MAX_LOOPS * step;
            log(threadId, "Building wheels");
            Wheel[] localWheels = initSieveWheels(position);
            log(threadId, "Started");

            Work ws = new Work(sieveVectorBound, estimateBucketCapacity());

            for (int loops = 0; loops < MAX_LOOPS; loops++) {
                position += step;

                sieveOnePass(position, localWheels, ws);

                // Merge batches
                if (ws.localBSmoothCount > 0) {
                    synchronized (bSmoothVectors) {
                        bSmoothVectors.addAll(ws.localBSmooth);
                    }
                    ws.localBSmooth.clear();
                    ws.localBSmoothCount = 0;
                }
                if (!ws.localBigPrimes.isEmpty()) {
                    synchronized (bigPrimesList) {
                        for (Pair p : ws.localBigPrimes) bigPrimesList.add(p.prime, p.data);
                    }
                    ws.localBigPrimes.clear();
                }

                // stats
                speed.incrementAndGet();
                if (speedCounter.incrementAndGet() == LOGS_TIME_BY_LOOPS) {
                    speedCounter.set(0);
                    logProcesses();
                }

                // solve?
                if (isReadyToBeSolved()) {
                    log(threadId, "Getting ready to solve");
                    if (threadId != 0) {
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
                        System.exit(0);
                        break;
                    } else {
                        synchronized (QuadraticThieve.this) {
                            QuadraticThieve.this.notifyAll();
                        }
                        System.exit(0);
                    }
                }
            }

            basePosition += (long) step * MAX_LOOPS * threadCount;
        }
    }

    private void logProcesses() {
        int loopsDone = this.speed.intValue();
        long now = System.currentTimeMillis();
        long sec = Math.max(1, (now - startingTime) / 1000);
        long blocksPerSec = loopsDone / sec;
        long valuesPerSec = blocksPerSec * (long) step;

        log("speed", valuesPerSec / 1000, "kValues a second, B-Smooth found",
                bSmoothVectors.size(), "Big primes found", bigPrimesList.getPrimesFound());
    }

    private Wheel[] initSieveWheels(long position) {
        Wheel[] wheels = new Wheel[B_SMOOTH];
        for (int i = 0; i < wheels.length; i++) {
            wheels[i] = new Wheel(primeBase[i], N, root.add(BigInteger.valueOf(position)), sieveVectorBound);
        }
        return wheels;
    }

    // Heuristic capacity: ~ block * (2 * sum_{p<=P} 1/p)
    private int estimateBucketCapacity() {
        // For B=5000 up to ~1e5, sum 1/p ~ 2.2–2.5; times 2 roots → ~4.5–5.
        // With block ~1e5 → ~5e5 hits. Add margin.
        return sieveVectorBound * 6;
    }

    /**
     * Single-pass sieve with buckets:
     *  - accumulate log(p) in acc[idx]
     *  - record prime hits into per-index linked lists
     *  - candidate scan by exact expected log
     *  - TD only against primes in that index's bucket
     */
    private void sieveOnePass(long destination, Wheel[] wheels, Work ws) {
        Arrays.fill(ws.acc, 0);
        Arrays.fill(ws.head, -1);
        ws.ptr = 0;

        // accumulate logs and fill buckets
        for (int i = 0; i < primeBaseInt.length; i++) {
            Wheel w = wheels[i];
            w.prepareToMove();
            int inc = logPScaled[i];
            while (w.testMove()) {
                int idx = w.move();
                if (idx >= 0 && idx < sieveVectorBound) {
                    ws.acc[idx] += inc;
                    // push (i) into bucket for idx
                    if (ws.ptr == ws.who.length) ws.growBuckets();
                    ws.who[ws.ptr] = i;
                    ws.next[ws.ptr] = ws.head[idx];
                    ws.head[idx] = ws.ptr;
                    ws.ptr++;
                }
            }
        }

        // scan candidates
        for (int idx = 0; idx < sieveVectorBound; idx++) {
            long tLong = destination + idx - sieveVectorBound;
            double expectedLn = Math.log(Math.max(1.0, (double) tLong * (tLong + double2Root)));
            int expectedScaled = (int) Math.round(expectedLn * LOG_SCALE);
            int remainderScaled = expectedScaled - ws.acc[idx];

            if (remainderScaled <= BIG_PRIME_CUTOFF_SCALED + LOG_EPS) {
                trialDivideBucketed(tLong, ws.head[idx], ws, expectedScaled);
            }
        }
    }

    /**
     * Divide Q(t) only by the primes that actually hit this index (bucket list).
     * Handles prime powers by repeated division; builds parity vector.
     * Accepts:
     *  - B-smooth (rem == 1)
     *  - single large prime (rem <= B_max^2 and prime)
     */
    private void trialDivideBucketed(long tLong, int headNode, Work ws, int expectedScaled) {
        if (headNode == -1) return;

        BigInteger t = BigInteger.valueOf(tLong);
        BigInteger x = root.add(t);
        BigInteger Q = x.multiply(x).subtract(N).abs();
        if (Q.signum() == 0) return;

        BitSet bits = new BitSet(B_SMOOTH);
        BigInteger rem = Q;

        // only the primes that hit this position
        for (int node = headNode; node != -1; node = ws.next[node]) {
            int i = ws.who[node];
            int p = primeBaseInt[i];

            if (p == 2) {
                int parity = 0;
                while (!rem.testBit(0)) {
                    rem = rem.shiftRight(1);
                    parity ^= 1;
                }
                if (parity == 1) bits.set(i);
            } else {
                int parity = 0;
                BigInteger bp = BigInteger.valueOf(p);
                // divide out p fully
                while (rem.mod(bp).signum() == 0) {
                    rem = rem.divide(bp);
                    parity ^= 1;
                }
                if (parity == 1) bits.set(i);
            }

            if (rem.equals(BigInteger.ONE)) break;
        }

        // tiny early abort: if rem is already huge, bail quickly
        if (!rem.equals(BigInteger.ONE) && rem.bitLength() - largePrimeBoundBI.bitLength() > 2) {
            return;
        }

        VectorData vd = new VectorData(bits, tLong);

        if (rem.equals(BigInteger.ONE)) {
            ws.localBSmooth.add(vd);
            ws.localBSmoothCount++;
            return;
        }

        if (rem.compareTo(largePrimeBoundBI) <= 0 && rem.isProbablePrime(20)) {
            ws.localBigPrimes.add(new Pair(rem.longValue(), vd));
        }
    }

    private boolean tryToSolve() {
        log("Building matrix");
        ArrayList<VectorData> vectorDatas;
        synchronized (bSmoothVectors) {
            vectorDatas = vectorsShrinker.shrink(bSmoothVectors, bigPrimesList);
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
        int i = 0;
        while (i < B_SMOOTH) {
            prime = prime.nextProbablePrime();
            if (MathUtils.isRootInQuadraticResidues(N, prime)) {
                primeBase[i] = prime;
                i++;
            }
        }
    }

    // --- thread-local work buffers ---
    private static final class Work {
        final int[] acc;              // scaled log sums
        final int[] head;             // per-index linked list head
        int[] who;                    // prime index for node
        int[] next;                   // next node pointer
        int ptr = 0;                  // current fill pointer

        final ArrayList<VectorData> localBSmooth = new ArrayList<>(256);
        final ArrayList<Pair> localBigPrimes = new ArrayList<>(256);
        int localBSmoothCount = 0;

        Work(int n, int bucketCap) {
            this.acc = new int[n];
            this.head = new int[n];
            this.who = new int[bucketCap];
            this.next = new int[bucketCap];
        }

        void growBuckets() {
            int newCap = who.length * BUCKET_GROWTH;
            who = Arrays.copyOf(who, newCap);
            next = Arrays.copyOf(next, newCap);
        }
    }

    private static final class Pair {
        final long prime;
        final VectorData data;
        Pair(long p, VectorData d) { this.prime = p; this.data = d; }
    }
}
