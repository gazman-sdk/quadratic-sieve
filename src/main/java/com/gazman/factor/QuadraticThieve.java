package com.gazman.factor;

import com.gazman.factor.matrix.BitMatrix;
import com.gazman.factor.wheels.Wheel;
import com.gazman.math.MathUtils;
import com.gazman.math.SqrRoot;

import java.math.BigInteger;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Quadratic Sieve core with:
 * - single-pass bucketed sieve
 * - k-way (k=5) large-prime hyper-merge (accept up to 4 LPs, merge aggressively)
 *
 * This replaces the previous 2-LP BigPrimePairs/BigPrimesList path.
 * Next step (separate patch): plug in the "search-aided pre-credit" residue scheduler
 * described in the paper to bias indices and shrink the remainder band.
 */
public class QuadraticThieve extends Logger {
    // ===== Tuning =====
    private static final int B_SMOOTH = 5000;               // factor base size
    public static final int MAX_LOOPS = B_SMOOTH * 2;       // blocks per thread before moving window
    public static final int LOGS_TIME_BY_LOOPS = B_SMOOTH / 20;

    private static final int SMALL_PRIME_LOG_ONLY_CUTOFF = 256;

    private static final int LOG_SCALE = 256;
    private static final int LOG_EPS = 6;

    private static final int BUCKET_GROWTH = 2;

    // k-split knob (we allow up to k-1 leftover large primes)
    private static final int K_SPLIT = 5;
    private static final int MAX_LP_COUNT = K_SPLIT - 1;

    // Accept only "small" leftovers for fast 64-bit factorization (sweet spot)
    private static final int MAX_RESIDUAL_BITS = 62; // < 2^62

    private final BigInteger N;
    private final BigInteger root;
    private final double double2Root;

    private final BigInteger[] primeBase = new BigInteger[B_SMOOTH];
    private final int[] primeBaseInt = new int[B_SMOOTH];
    private final int[] logPScaled = new int[B_SMOOTH];

    private final int sieveVectorBound; // block size
    private final int step;

    private final int threadCount = 2;

    private final ArrayList<VectorData> bSmoothVectors = new ArrayList<>();

    private final AtomicInteger speedCounter = new AtomicInteger(0);
    private final AtomicInteger speed = new AtomicInteger(0);

    private final int LN_2ROOT_SCALED;

    private volatile long startingTime;

    // New: k-way hyper-merge for up to 4 LPs
    private final LargePrimeKMerger kMerger;

    public QuadraticThieve(BigInteger input) {
        log("Factoring started");
        N = input;
        root = SqrRoot.bigIntSqRootCeil(input);
        double2Root = root.add(root).doubleValue();
        LN_2ROOT_SCALED = (int) Math.round(Math.log(double2Root) * LOG_SCALE);

        log("Building Prime Base");
        buildPrimeBase();

        BigInteger highestPrime = primeBase[primeBase.length - 1];
        sieveVectorBound = highestPrime.intValue();
        step = sieveVectorBound;

        for (int i = 0; i < B_SMOOTH; i++) {
            int p = primeBase[i].intValue();
            primeBaseInt[i] = p;
            logPScaled[i] = (int) Math.round(Math.log(p) * LOG_SCALE);
        }

        kMerger = new LargePrimeKMerger(N);

        log("Biggest prime is", highestPrime);
        log();
        log("Working on", threadCount, "threads");
        log("Start searching");
    }

    private static int fastLogScaled(long x) {
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

                // Merge thread-local batches into global state
                if (ws.localBSmoothCount > 0) {
                    synchronized (bSmoothVectors) {
                        bSmoothVectors.addAll(ws.localBSmooth);
                    }
                    ws.localBSmooth.clear();
                    ws.localBSmoothCount = 0;
                }
                if (!ws.localPartials.isEmpty()) {
                    for (Partial p : ws.localPartials) {
                        VectorData merged = kMerger.add(p.primes, p.data);
                        if (merged != null) {
                            synchronized (bSmoothVectors) {
                                bSmoothVectors.add(merged);
                            }
                        }
                    }
                    ws.localPartials.clear();
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
                bSmoothVectors.size());
    }

    private Wheel[] initSieveWheels(long position) {
        Wheel[] wheels = new Wheel[B_SMOOTH];
        for (int i = 0; i < wheels.length; i++) {
            wheels[i] = new Wheel(primeBase[i], N, root.add(BigInteger.valueOf(position)), sieveVectorBound);
        }
        return wheels;
    }

    private int estimateBucketCapacity() {
        return sieveVectorBound * 3;
    }

    /**
     * Single-pass bucketed sieve as before; TD changed to support k-way LP handling.
     */
    private void sieveOnePass(long destination, Wheel[] wheels, Work ws) {
        Arrays.fill(ws.acc, 0);
        Arrays.fill(ws.head, -1);
        ws.ptr = 0;

        // Pass 1: tiny primes (log only)
        for (int i = 0; i < primeBaseInt.length; i++) {
            if (primeBaseInt[i] >= SMALL_PRIME_LOG_ONLY_CUTOFF) break;
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

        // Pass 2: larger primes (logs + per-index bucket lists)
        for (int i = 0; i < primeBaseInt.length; i++) {
            if (primeBaseInt[i] < SMALL_PRIME_LOG_ONLY_CUTOFF) continue;
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

        long t0 = Math.max(1L, destination - sieveVectorBound);
        int expectedBlockScaled = LN_2ROOT_SCALED + fastLogScaled(t0);

        final int cutoff = (int) Math.round(Math.log(Math.pow(primeBase[primeBase.length - 1].doubleValue(), 2.0)) * LOG_SCALE) + LOG_EPS;
        for (int idx = 0; idx < sieveVectorBound; idx++) {
            int remainderScaled = expectedBlockScaled - ws.acc[idx];
            if (remainderScaled <= cutoff) {
                long tLong = destination + idx - sieveVectorBound;
                trialDivideBucketed(tLong, ws.head[idx], ws);
            }
        }
    }

    /**
     * Trial divide Q(t) = (root+t)^2 - N by:
     *  - tiny primes (pass 1)
     *  - bucketed larger primes (pass 2 list)
     *
     * Accept:
     *  - B-smooth (remainder 1 or perfect square)
     *  - residue that fits in 62 bits and is product of <= 4 primes → k-way merge
     */
    private void trialDivideBucketed(long tLong, int headNode, Work ws) {
        BigInteger t = BigInteger.valueOf(tLong);
        BigInteger x = root.add(t);
        BigInteger Q = x.multiply(x).subtract(N).abs();
        if (Q.signum() == 0) return;

        BitSet bits = new BitSet(B_SMOOTH);
        BigInteger rem = Q;

        // Divide by small primes (not bucketed)
        for (int i = 0; i < primeBaseInt.length; i++) {
            int p = primeBaseInt[i];
            if (p >= SMALL_PRIME_LOG_ONLY_CUTOFF) break;
            if (p == 2) {
                int twoFactors = rem.getLowestSetBit();
                if (twoFactors > 0) {
                    rem = rem.shiftRight(twoFactors);
                    if ((twoFactors & 1) == 1) bits.set(i);
                }
            } else {
                int parity = 0;
                BigInteger pBigInt = primeBase[i];
                while (true) {
                    BigInteger[] dr = rem.divideAndRemainder(pBigInt);
                    if (dr[1].signum() == 0) {
                        rem = dr[0];
                        parity ^= 1;
                    } else break;
                }
                if (parity == 1) bits.set(i);
            }
        }

        // Divide by the larger primes from the bucket list.
        for (int node = headNode; node != -1; node = ws.next[node]) {
            int i = ws.who[node];
            int parity = 0;
            BigInteger pBigInt = primeBase[i];
            while (true) {
                BigInteger[] dr = rem.divideAndRemainder(pBigInt);
                if (dr[1].signum() == 0) {
                    rem = dr[0];
                    parity ^= 1;
                } else break;
            }
            if (parity == 1) bits.set(i);
        }

        // Fully B-smooth
        if (rem.equals(BigInteger.ONE)) {
            VectorData vd = new VectorData(bits, tLong);
            vd.x = x;      // keep x, y for congruence
            vd.y = Q;
            ws.localBSmooth.add(vd);
            ws.localBSmoothCount++;
            return;
        }

        // If the remainder is a perfect square > 1, it's also OK (parity-wise)
        BigInteger sq = SqrRoot.bigIntSqRootFloor(rem);
        if (sq.multiply(sq).equals(rem)) {
            VectorData vd = new VectorData(bits, tLong);
            vd.x = x;
            vd.y = Q;
            ws.localBSmooth.add(vd);
            ws.localBSmoothCount++;
            return;
        }

        // Else: try k-way LP handling if remainder fits into 62 bits
        if (rem.bitLength() > MAX_RESIDUAL_BITS) return;

        long r = rem.longValue(); // safe (by bitLength check)
        long[] oddPrimes = factorOddResidueUpTo4(r);
        if (oddPrimes == null || oddPrimes.length == 0 || oddPrimes.length > MAX_LP_COUNT) return;

        VectorData cur = new VectorData(bits, tLong);
        cur.x = x;
        cur.y = Q;
        ws.localPartials.add(new Partial(oddPrimes, cur));
    }

    private boolean tryToSolve() {
        log("Building matrix");
        ArrayList<VectorData> vectorDatas;
        synchronized (bSmoothVectors) {
            vectorDatas = new ArrayList<>(bSmoothVectors);
        }

        BitMatrix bitMatrix = new BitMatrix();
        ArrayList<ArrayList<VectorData>> solutions = bitMatrix.solve(vectorDatas);

        for (int i = 0; i < solutions.size(); i++) {
            ArrayList<VectorData> solution = solutions.get(i);
            if (solution.isEmpty()) continue;
            log("Testing solution", (i + 1) + "/" + solutions.size());
            if (testSolution(solution)) return true;
        }
        log("no luck");
        return false;
    }

    private boolean isReadyToBeSolved() {
        // with hyper-merge we wait for a little margin beyond base size
        return bSmoothVectors.size() >= B_SMOOTH + 10;
    }

    private boolean testSolution(ArrayList<VectorData> solutionVector) {
        BigInteger one = BigInteger.ONE;
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

    // ===== Thread-local work buffers =====
    private static final class Work {
        final int[] acc;
        final int[] head;
        final ArrayList<VectorData> localBSmooth = new ArrayList<>(256);
        final ArrayList<Partial> localPartials = new ArrayList<>(256);
        int[] who;
        int[] next;
        int ptr = 0;
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

    private record Partial(long[] primes, VectorData data) {}

    // ===== 64-bit factoring (fast path for small remnants) =====

    /**
     * Factor a 62-bit number and return the DISTINCT primes with ODD exponent, sorted.
     * Return null if factoring fails quickly or if exponent parities cancel entirely (handled elsewhere).
     */
    private static long[] factorOddResidueUpTo4(long n) {
        if (n <= 1) return new long[0];
        Map<Long, Integer> mp = new HashMap<>();
        factor64(n, mp);

        // Collect primes with odd exponent
        ArrayList<Long> list = new ArrayList<>(4);
        for (Map.Entry<Long, Integer> e : mp.entrySet()) {
            if ((e.getValue() & 1) == 1) list.add(e.getKey());
        }
        if (list.isEmpty()) return new long[0];
        if (list.size() > 4) return null;
        list.sort(Comparator.naturalOrder());
        long[] out = new long[list.size()];
        for (int i = 0; i < out.length; i++) out[i] = list.get(i);
        return out;
    }

    private static void factor64(long n, Map<Long, Integer> acc) {
        if (n <= 1) return;
        if ((n & 1L) == 0L) {
            int c = 0; while ((n & 1L) == 0L) { n >>>= 1; c++; }
            acc.merge(2L, c, Integer::sum);
            if (n > 1) factor64(n, acc);
            return;
        }
        if (isPrime64(n)) {
            acc.merge(n, 1, Integer::sum);
            return;
        }
        long d = pollardRho64(n);
        if (d == n) { // fallback (rare)
            // try small trial division up to 100k
            for (long p = 3; p <= 100_000 && p * p <= n; p += 2) {
                if (n % p == 0) {
                    int c = 0; while (n % p == 0) { n /= p; c++; }
                    acc.merge(p, c, Integer::sum);
                    factor64(n, acc);
                    return;
                }
            }
            // give up: treat n as prime (should be very rare with 62-bit cap)
            acc.merge(n, 1, Integer::sum);
            return;
        }
        factor64(d, acc);
        factor64(n / d, acc);
    }

    // Deterministic MR for 64-bit
    private static boolean isPrime64(long n) {
        if (n < 2) return false;
        for (long p : new long[]{2,3,5,7,11,13,17,19,23,29,31,37}) {
            if (n % p == 0) return n == p;
        }
        long d = n - 1, s = 0;
        while ((d & 1L) == 0L) { d >>>= 1; s++; }
        long[] bases = {2, 3, 5, 7, 11, 13};
        for (long a : bases) {
            if (a % n == 0) continue;
            long x = powMod64(a, d, n);
            if (x == 1 || x == n - 1) continue;
            boolean cont = false;
            for (int r = 1; r < s; r++) {
                x = mulMod64(x, x, n);
                if (x == n - 1) { cont = true; break; }
            }
            if (!cont) return false;
        }
        return true;
    }

    private static long pollardRho64(long n) {
        if ((n & 1L) == 0L) return 2L;
        java.util.Random rnd = new java.util.Random(42);
        while (true) {
            long c = 1 + Math.abs(rnd.nextLong()) % (n - 1);
            long x = 2 + Math.abs(rnd.nextLong()) % (n - 2);
            long y = x;
            long d = 1;
            while (d == 1) {
                x = (mulMod64(x, x, n) + c) % n;
                y = (mulMod64(y, y, n) + c) % n;
                y = (mulMod64(y, y, n) + c) % n;
                long diff = x > y ? x - y : y - x;
                d = gcd64(diff, n);
            }
            if (d != n) return d;
        }
    }

    private static long gcd64(long a, long b) {
        while (b != 0) {
            long t = a % b;
            a = b; b = t;
        }
        return a;
    }

    // 128-bit safe via BigInteger (fast enough for 62-bit inputs)
    private static long mulMod64(long a, long b, long mod) {
        return BigInteger.valueOf(a).multiply(BigInteger.valueOf(b)).mod(BigInteger.valueOf(mod)).longValue();
    }
    private static long powMod64(long a, long e, long mod) {
        return BigInteger.valueOf(a).modPow(BigInteger.valueOf(e), BigInteger.valueOf(mod)).longValue();
    }
}
