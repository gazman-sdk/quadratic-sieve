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
 * - Single-pass sieve (scaled logs)
 * - Per-index bucket lists of hitting primes (massive TD reduction)
 * - Single-large-prime policy aligned to original (<= B_max^2)
 * - Per-thread batching (no hot-path sync/alloc)
 * - Cache-optimized hybrid sieving for small/large primes
 */
public class QuadraticThieve extends Logger {
    // --- Tuning knobs ---
    private static final int B_SMOOTH = 5000;               // factor base size
    public static final int MAX_LOOPS = B_SMOOTH * 2;       // blocks per thread before moving window
    public static final int LOGS_TIME_BY_LOOPS = B_SMOOTH / 20;

    // Primes smaller than this are sieved without being added to the expensive
    // per-index bucket lists, drastically improving cache performance.
    private static final int SMALL_PRIME_LOG_ONLY_CUTOFF = 256;

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

    private final int cores = Runtime.getRuntime().availableProcessors();
    // Threads
    private final int threadCount = cores;

    private final VectorsShrinker vectorsShrinker = new VectorsShrinker();
    private final BigPrimesList bigPrimesList = new BigPrimesList();

    private final ArrayList<VectorData> bSmoothVectors = new ArrayList<>();

    private final AtomicInteger speedCounter = new AtomicInteger(0);
    private final AtomicInteger speed = new AtomicInteger(0);
    // --- fast log threshold ---
    private final int LN_2ROOT_SCALED;          // round( ln(2*sqrt(N)) * LOG_SCALE )
    // 2-large-prime aggregator
    private final BigPrimePairs bigPrimePairs;
    private volatile long startingTime;


    public QuadraticThieve(BigInteger input) {
        log("Factoring started");
        N = input;
        root = SqrRoot.bigIntSqRootCeil(input);
        double2Root = root.add(root).doubleValue();
        LN_2ROOT_SCALED = (int) Math.round(Math.log(double2Root) * LOG_SCALE);
        bigPrimePairs = new BigPrimePairs(root, N);

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

    private static int fastLogScaled(long x) {
        // one Math.log per block; here keep it inline for clarity
        return (int) Math.round(Math.log((double) x) * LOG_SCALE);
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
        // We only bucket large primes, so capacity can be smaller.
        // Heuristic: ~ block * (2 * sum_{p > cutoff} 1/p)
        return sieveVectorBound * 3;
    }

    /**
     * Single-pass sieve with buckets:
     * - accumulate log(p) in acc[idx]
     * - record prime hits into per-index linked lists
     * - candidate scan by exact expected log
     * - TD only against primes in that index's bucket
     */
    private void sieveOnePass(long destination, Wheel[] wheels, Work ws) {
        Arrays.fill(ws.acc, 0);
        Arrays.fill(ws.head, -1);
        ws.ptr = 0;

        // --- OPTIMIZATION ---
        // Pass 1: Sieve small primes (logs only) for cache performance.
        for (int i = 0; i < primeBaseInt.length; i++) {
            if (primeBaseInt[i] >= SMALL_PRIME_LOG_ONLY_CUTOFF) {
                break;
            }
            Wheel w = wheels[i];
            w.prepareToMove();
            int inc = logPScaled[i];
            while (w.testMove()) {
                int idx = w.move();
                if (idx >= 0 && idx < sieveVectorBound) {
                    ws.acc[idx] += inc;
                }
            }
        }

        // Pass 2: Sieve larger primes (logs and bucket lists).
        // This is much faster as these primes have far fewer hits.
        for (int i = 0; i < primeBaseInt.length; i++) {
            if (primeBaseInt[i] < SMALL_PRIME_LOG_ONLY_CUTOFF) {
                continue;
            }
            Wheel w = wheels[i];
            w.prepareToMove();
            int inc = logPScaled[i];
            while (w.testMove()) {
                int idx = w.move();
                if (idx >= 0 && idx < sieveVectorBound) {
                    ws.acc[idx] += inc;
                    if (ws.ptr == ws.who.length) ws.growBuckets();
                    ws.who[ws.ptr] = i;
                    ws.next[ws.ptr] = ws.head[idx];
                    ws.head[idx] = ws.ptr;
                    ws.ptr++;
                }
            }
        }


        // --- block-constant expected log (no per-index Math.log) ---
        // t runs from (destination - sieveVectorBound) ... (destination - 1)
        long t0 = Math.max(1L, destination - sieveVectorBound); // avoid log(0)
        int expectedBlockScaled = LN_2ROOT_SCALED + fastLogScaled(t0);

        // scan candidates using a *constant* threshold for this block
        final int cutoff = BIG_PRIME_CUTOFF_SCALED + LOG_EPS;
        for (int idx = 0; idx < sieveVectorBound; idx++) {
            int remainderScaled = expectedBlockScaled - ws.acc[idx];
            if (remainderScaled <= cutoff) {
                long tLong = destination + idx - sieveVectorBound;
                trialDivideBucketed(tLong, ws.head[idx], ws);
            }
        }
    }

    /**
     * Divide Q(t) only by the primes that actually hit this index (bucket list).
     * Handles prime powers by repeated division; builds parity vector.
     * Accepts:
     * - B-smooth (rem == 1)
     * - single large prime (rem <= B_max^2 and prime)
     */
    private void trialDivideBucketed(long tLong, int headNode, Work ws) {
        BigInteger t = BigInteger.valueOf(tLong);
        BigInteger x = root.add(t);
        BigInteger Q = x.multiply(x).subtract(N).abs();
        if (Q.signum() == 0) return;

        BitSet bits = new BitSet(B_SMOOTH);
        BigInteger rem = Q;

        // Step 1: Trial divide by the small primes that were not bucketed.
        for (int i = 0; i < primeBaseInt.length; i++) {
            int p = primeBaseInt[i];
            if (p >= SMALL_PRIME_LOG_ONLY_CUTOFF) {
                break;
            }

            if (p == 2) {
                int twoFactors = rem.getLowestSetBit();
                if (twoFactors > 0) {
                    rem = rem.shiftRight(twoFactors);
                    if ((twoFactors & 1) == 1) {
                        bits.set(i);
                    }
                }
            } else {
                int parity = 0;
                BigInteger pBigInt = primeBase[i];
                while (true) {
                    BigInteger[] divRem = rem.divideAndRemainder(pBigInt);
                    if (divRem[1].signum() == 0) {
                        rem = divRem[0];
                        parity ^= 1;
                    } else {
                        break;
                    }
                }
                if (parity == 1) bits.set(i);
            }
        }


        // Step 2: Divide by the larger primes from the bucket list.
        for (int node = headNode; node != -1; node = ws.next[node]) {
            int i = ws.who[node];
            int parity = 0;
            BigInteger pBigInt = primeBase[i];
            while (true) {
                BigInteger[] divRem = rem.divideAndRemainder(pBigInt);
                if (divRem[1].signum() == 0) {
                    rem = divRem[0];
                    parity ^= 1;
                } else {
                    break;
                }
            }
            if (parity == 1) bits.set(i);
        }

        // fully B-smooth
        if (rem.equals(BigInteger.ONE)) {
            ws.localBSmooth.add(new VectorData(bits, tLong));
            ws.localBSmoothCount++;
            return;
        }

        // Only consider small leftovers (<= B_max^2)
        if (rem.compareTo(largePrimeBoundBI) > 0) return;

        // --- Single large prime path (existing behavior) ---
        if (rem.isProbablePrime(20)) {
            ws.localBigPrimes.add(new Pair(rem.longValue(), new VectorData(bits, tLong)));
            return;
        }

        // --- 2-large-prime path: try to split rem into p*q with p,q <= B_max^2 ---
        long r = rem.longValue(); // safe: rem <= B_max^2 ~ 1e10
        PairLong pq = factorSemiprimeLE1e10(r);
        if (pq == null) return; // not of the right shape; ignore

        // both factors prime?
        if (!BigInteger.valueOf(pq.a).isProbablePrime(20)) return;
        if (!BigInteger.valueOf(pq.b).isProbablePrime(20)) return;

        // Try to merge with a matching (p,q) relation to produce B-smooth immediately
        VectorData cur = new VectorData(bits, tLong);
        VectorData merged = bigPrimePairs.add(pq.a, pq.b, cur);
        if (merged != null) {
            ws.localBSmooth.add(merged);
            ws.localBSmoothCount++;
        }
    }

    // rem <= B_max^2 ~ 1e10 : at least one factor <= 1e5
    private PairLong factorSemiprimeLE1e10(long n) {
        if ((n & 1L) == 0) return new PairLong(2L, n >>> 1);
        // quick checks by 3 and 5
        if (n % 3L == 0) return new PairLong(3L, n / 3L);
        if (n % 5L == 0) return new PairLong(5L, n / 5L);

        // odd trial division up to 100,000
        for (long d = 7; d * d <= n && d <= 100_000L; d += 2) {
            if (n % d == 0) {
                return new PairLong(d, n / d);
            }
        }
        return null;
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
            if (solution.isEmpty()) {
                continue;
            }
            log("Testing solution", (i + 1) + "/" + solutions.size());
            if (testSolution(solution)) {
                return true;
            }
        }
        log("no luck");
        return false;
    }

    private boolean isReadyToBeSolved() {
        return bSmoothVectors.size() + bigPrimesList.getPrimesFound() >= B_SMOOTH + 10; // Add a margin
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
        BigInteger gcd = N.gcd(x.subtract(y).abs());
        if (!gcd.equals(one) && !gcd.equals(N)) {
            log("Solved");
            log(gcd);
            return true;
        }

        gcd = N.gcd(x.add(y));
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
        final ArrayList<VectorData> localBSmooth = new ArrayList<>(256);
        final ArrayList<Pair> localBigPrimes = new ArrayList<>(256);
        int[] who;                    // prime index for node
        int[] next;                   // next node pointer
        int ptr = 0;                  // current fill pointer
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

    private record PairLong(long a, long b) {
    }

    private record Pair(long prime, VectorData data) {
    }
}