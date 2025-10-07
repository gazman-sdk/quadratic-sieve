package com.gazman.factor;

import com.gazman.factor.matrix.BitMatrix;
import com.gazman.factor.matrix.VectorsShrinker;
import com.gazman.math.MathUtils;
import com.gazman.math.SqrRoot;
import com.gazman.search.FSTSearch;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * QS refactored to k=5 linear product with FSTSearch.
 *
 * CHANGE IN THIS REV:
 *  - Independent sieving of the 4 moving linears (L1..L4) instead of screening on the product log.
 *  - Four per-form accumulators, four per-form residue wheels, and per-block intersection stats.
 *  - TD runs only when ALL FOUR forms pass their own log threshold at the same index.
 *
 * Forms (moving part only):
 *   L1(X) = |sX - t|
 *   L2(X) = |sX + t|
 *   L3(X) = |f s^2 X - c|
 *   L4(X) = |f s^2 X + c|
 * where X = M + x  (M is search center; x is sieve offset)
 */
public class QuadraticThieve extends Logger {

    // ===== Debug / logging knobs =====
    private static final boolean DEBUG_VERBOSE = true;    // richer per-block logs
    private static final int     DEBUG_BLOCKS  = 3;       // log first N blocks verbosely on thread 0
    private static final long    DEBUG_EVERY_MS = 5000;   // then every N ms

    private static final int     SHOW_TOP = 8;            // leaderboard size (product-style, retained)
    private static final boolean REQUIRE_ALL_FORMS = true;// gate TD on all four passing form thresholds

    // ===== Tuning =====
    private static final int B_SMOOTH = 5000;
    public  static final int MAX_LOOPS = B_SMOOTH * 2;
    public  static final int LOGS_TIME_BY_LOOPS = B_SMOOTH / 20;

    private static final int SMALL_PRIME_LOG_ONLY_CUTOFF = 256;

    private static final int LOG_SCALE = 256;  // natural logs scaled
    private static final int LOG_EPS   = 6;
    private static final int BUCKET_GROWTH = 2;

    // ===== Problem / factor base =====
    private final BigInteger N;

    private final BigInteger[] primeBase = new BigInteger[B_SMOOTH];
    private final int[]        primeBaseInt = new int[B_SMOOTH];
    private final int[]        logPScaled   = new int[B_SMOOTH];

    private final int sieveVectorBound; // block size
    private final int step;

    private final BigInteger largePrimeBoundBI;
    private final int BIG_PRIME_CUTOFF_SCALED;

    private final int threadCount = Math.max(1, Runtime.getRuntime().availableProcessors());

    private final VectorsShrinker vectorsShrinker = new VectorsShrinker();
    private final BigPrimesList  bigPrimesList    = new BigPrimesList();
    private final ArrayList<VectorData> bSmoothVectors = new ArrayList<>();

    private final AtomicInteger speedCounter = new AtomicInteger(0);
    private final AtomicInteger speed        = new AtomicInteger(0);

    // ===== k=5 polynomial from FSTSearch =====
    private final FSTSearch.Poly poly;
    private final BigInteger f, s, t, c;
    private final long       M;      // search center
    private final BigInteger MB;     // as BigInteger

    private volatile long startingTime;
    private volatile long lastDebug = 0;

    // Pairing for 2-LP (product level)
    private final BigPrimePairs bigPrimePairs;

    private final Random rnd = new Random(42);

    public QuadraticThieve(BigInteger input) {
        log("Factoring started");
        this.N = input;

        log("Building factor base...");
        buildPrimeBase();

        BigInteger highestPrime = primeBase[B_SMOOTH - 1];
        this.sieveVectorBound = highestPrime.intValue();
        this.step = sieveVectorBound;

        for (int i = 0; i < B_SMOOTH; i++) {
            int p = primeBase[i].intValue();
            primeBaseInt[i] = p;
            logPScaled[i] = (int) Math.round(Math.log(p) * LOG_SCALE);
        }

        // Choose M ~ 2^(|N|/5), capped at 2^60 to stay in long
        int portionBits = Math.max(1, N.bitLength() / 5);
        int mBits = Math.min(60, portionBits);
        long searchM = 1L << mBits;
        this.M  = searchM;
        this.MB = BigInteger.valueOf(M);

        // Use FSTSearch with time you set inside FSTSearch (you asked to keep it constant there)
        FSTSearch search = new FSTSearch(N);
        this.poly = search.findBestPolynomial(searchM);
        this.f = poly.f; this.s = poly.s; this.t = poly.t; this.c = poly.c;

        // thresholds on primes (LP bound)
        this.largePrimeBoundBI = highestPrime.multiply(highestPrime);
        this.BIG_PRIME_CUTOFF_SCALED =
                (int) Math.round(Math.log(largePrimeBoundBI.doubleValue()) * LOG_SCALE);

        // init shrinker and LP merger for k=5 semantics (same as previous rev)
        vectorsShrinker.initK5(primeBase.length, N, f, s, t, c, MB);
        bigPrimePairs = new BigPrimePairs(N, f, s, t, c, MB);

        log("k=5 polynomial:");
        log("  f =", f, "s =", s, "t =", t, "c =", c, " | score:", String.format("%.2f", poly.score));
        log("M =", M, " (targetBits≈", portionBits, ")");
        log("Factor base:", B_SMOOTH, " | biggest prime:", highestPrime);
        logCenteringFeasibility();  // tells if this poly can in-principle be centered
        calibrateAndLogOnce();
        log();
        log("Working on", threadCount, "threads");
        log("Start sieving (k=5; independent forms)...");
    }

    public void start() {
        for (int i = 0; i < threadCount; i++) {
            final int threadId = i;
            new Thread(() -> execute(threadId), "QS-" + threadId).start();
        }
    }

    private void execute(int threadId) {
        long basePosition = 0;
        if (threadId == 0) {
            startingTime = System.currentTimeMillis();
            lastDebug = startingTime;
        }

        int blocksSeen = 0;

        while (true) {
            long position = basePosition + (long) threadId * MAX_LOOPS * step;
            if (blocksSeen == 0) log(threadId, "Preparing linear wheels (k=5; per-form)");

            // Build 4 wheels per prime (L1..L4); M baked in
            SingleWheel[][] formWheels = initFormWheels();

            if (blocksSeen == 0) {
                wheelSanity(formWheels); // small sanity check on tiny primes
            }

            log(threadId, "Started");

            Work ws = new Work(sieveVectorBound, estimateBucketCapacity());

            for (int loops = 0; loops < MAX_LOOPS; loops++) {
                position += step;
                blocksSeen++;

                sieveOnePass(position, formWheels, ws, threadId, blocksSeen);

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

                // Solve?
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

        log("speed", valuesPerSec / 1000, "kValues/s, B-smooth", bSmoothVectors.size(),
                "Big primes", bigPrimesList.getPrimesFound());
    }

    // === NEW: 4 wheels per prime ===
    private SingleWheel[][] initFormWheels() {
        SingleWheel[][] wheels = new SingleWheel[B_SMOOTH][4];
        for (int i = 0; i < B_SMOOTH; i++) {
            int p = primeBaseInt[i];
            wheels[i][0] = SingleWheel.build(p, s,  t, MB);            // L1: sX ≡  t (mod p)
            wheels[i][1] = SingleWheel.build(p, s,  t.negate(), MB);   // L2: sX ≡ -t
            BigInteger fs2 = f.multiply(s).multiply(s);
            wheels[i][2] = SingleWheel.build(p, fs2, c, MB);           // L3: fs^2 X ≡  c
            wheels[i][3] = SingleWheel.build(p, fs2, c.negate(), MB);  // L4: fs^2 X ≡ -c
        }
        return wheels;
    }

    private int estimateBucketCapacity() { return sieveVectorBound * 3; }

    // === Single block pass with independent form accumulators ===
    private void sieveOnePass(long destination, SingleWheel[][] wheels, Work ws, int threadId, int blocksSeen) {
        Arrays.fill(ws.accL1, 0);
        Arrays.fill(ws.accL2, 0);
        Arrays.fill(ws.accL3, 0);
        Arrays.fill(ws.accL4, 0);

        Arrays.fill(ws.head, -1);
        ws.ptr = 0;

        long baseX = destination - sieveVectorBound; // x offset block start

        // Small primes: per-form logs only
        long smallHits = 0;
        for (int i = 0; i < primeBaseInt.length; i++) {
            if (primeBaseInt[i] >= SMALL_PRIME_LOG_ONLY_CUTOFF) break;
            int inc = logPScaled[i];

            for (int fIdx = 0; fIdx < 4; fIdx++) {
                SingleWheel w = wheels[i][fIdx];
                if (w == null) continue;
                w.rebase(baseX, sieveVectorBound);
                while (w.hasNext()) {
                    int idx = w.next();
                    if (idx < 0 || idx >= sieveVectorBound) continue;
                    switch (fIdx) {
                        case 0 -> ws.accL1[idx] += inc;
                        case 1 -> ws.accL2[idx] += inc;
                        case 2 -> ws.accL3[idx] += inc;
                        case 3 -> ws.accL4[idx] += inc;
                    }
                    smallHits++;
                }
            }
        }

        // Large primes: per-form logs + product-level buckets
        long bucketNodes = 0;
        for (int i = 0; i < primeBaseInt.length; i++) {
            if (primeBaseInt[i] < SMALL_PRIME_LOG_ONLY_CUTOFF) continue;
            int inc = logPScaled[i];

            for (int fIdx = 0; fIdx < 4; fIdx++) {
                SingleWheel w = wheels[i][fIdx];
                if (w == null) continue;
                w.rebase(baseX, sieveVectorBound);
                while (w.hasNext()) {
                    int idx = w.next();
                    if (idx < 0 || idx >= sieveVectorBound) continue;
                    switch (fIdx) {
                        case 0 -> ws.accL1[idx] += inc;
                        case 1 -> ws.accL2[idx] += inc;
                        case 2 -> ws.accL3[idx] += inc;
                        case 3 -> ws.accL4[idx] += inc;
                    }
                    if (ws.ptr == ws.who.length) ws.growBuckets();
                    ws.who[ws.ptr] = i;
                    ws.next[ws.ptr] = ws.head[idx];
                    ws.head[idx] = ws.ptr;
                    ws.ptr++;
                    bucketNodes++;
                }
            }
        }

        // Block-constant expected logs per form at center
        int[] expFormScaled = expectedFormLogsScaled(baseX + (sieveVectorBound >> 1));
        final int cutoff = BIG_PRIME_CUTOFF_SCALED + LOG_EPS; // ln*256
        final int sumCutoff = cutoff; // sum gate uses same bound (tune if desired)

        // Per-form candidate counts and intersections
        int cand1 = 0, cand2 = 0, cand3 = 0, cand4 = 0;
        int cand12 = 0, cand13 = 0, cand14 = 0, cand23 = 0, cand24 = 0, cand34 = 0;
        int cand123 = 0, cand124 = 0, cand134 = 0, cand234 = 0;
        int cand1234 = 0;
        int candSUM = 0; // NEW: sum-of-remainders ≤ cutoff

        // Legacy leaderboard (product style)
        int[] bestIdx = new int[SHOW_TOP];
        int[] bestRem = new int[SHOW_TOP];
        Arrays.fill(bestIdx, -1); Arrays.fill(bestRem, Integer.MAX_VALUE);
        int minRem = Integer.MAX_VALUE, maxRem = Integer.MIN_VALUE;

        int candidatesForTD = 0;

        // Scan indices
        for (int idx = 0; idx < sieveVectorBound; idx++) {
            int r1 = expFormScaled[0] - ws.accL1[idx];
            int r2 = expFormScaled[1] - ws.accL2[idx];
            int r3 = expFormScaled[2] - ws.accL3[idx];
            int r4 = expFormScaled[3] - ws.accL4[idx];

            boolean p1 = (r1 <= cutoff);
            boolean p2 = (r2 <= cutoff);
            boolean p3 = (r3 <= cutoff);
            boolean p4 = (r4 <= cutoff);

            if (p1) cand1++; if (p2) cand2++; if (p3) cand3++; if (p4) cand4++;
            if (p1 && p2) cand12++;
            if (p1 && p3) cand13++;
            if (p1 && p4) cand14++;
            if (p2 && p3) cand23++;
            if (p2 && p4) cand24++;
            if (p3 && p4) cand34++;
            if (p1 && p2 && p3) cand123++;
            if (p1 && p2 && p4) cand124++;
            if (p1 && p3 && p4) cand134++;
            if (p2 && p3 && p4) cand234++;
            if (p1 && p2 && p3 && p4) cand1234++;

            // NEW: sum-of-remainders gate
            int remSum = r1 + r2 + r3 + r4;
            if (remSum <= sumCutoff) {
                candSUM++;
                long x = destination + idx - sieveVectorBound; // x offset
                candidatesForTD++;
                trialDivideBucketed(x, ws.head[idx], ws);
            }

            // product-style tracking for intuition
            if (remSum < minRem) minRem = remSum;
            if (remSum > maxRem) maxRem = remSum;

            if (remSum < bestRem[SHOW_TOP - 1]) {
                int j = SHOW_TOP - 1;
                while (j > 0 && remSum < bestRem[j - 1]) {
                    bestRem[j] = bestRem[j - 1];
                    bestIdx[j] = bestIdx[j - 1];
                    j--;
                }
                bestRem[j] = remSum;
                bestIdx[j] = idx;
            }
        }

        // Emit diagnostics (first few blocks + periodic)
        boolean shouldLog = false;
        if (blocksSeen <= DEBUG_BLOCKS && threadId == 0) shouldLog = true;
        long now = System.currentTimeMillis();
        if (!shouldLog && (now - lastDebug) >= DEBUG_EVERY_MS && threadId == 0) {
            lastDebug = now; shouldLog = true;
        }

        if (shouldLog) {
            final double toBits = 1.0 / (Math.log(2.0) * LOG_SCALE);

            // per-form maxima over block
            long xStart = baseX, xEnd = baseX + sieveVectorBound - 1;
            int[] maxFormBits = formMaxBitsForBlock(xStart, xEnd);

            // “best-sum” index detail
            int bestI = bestIdx[0];
            long bestXoff = (bestI >= 0) ? (destination + bestI - sieveVectorBound) : (baseX + (sieveVectorBound >> 1));
            BigInteger Xbest = MB.add(BigInteger.valueOf(bestXoff));
            int[] bestFormBits = formBitsAtX(Xbest);
            int bestSumBits = bestFormBits[0] + bestFormBits[1] + bestFormBits[2] + bestFormBits[3];
            double bestAccBits =
                    (ws.accL1[(bestI >= 0) ? bestI : (sieveVectorBound >> 1)]
                            + ws.accL2[(bestI >= 0) ? bestI : (sieveVectorBound >> 1)]
                            + ws.accL3[(bestI >= 0) ? bestI : (sieveVectorBound >> 1)]
                            + ws.accL4[(bestI >= 0) ? bestI : (sieveVectorBound >> 1)]) * toBits;
            double bestRemBits = bestSumBits - bestAccBits;

            // Independence-based expectation for ALL4 (for intuition only)
            double p1f = cand1 / (double) sieveVectorBound;
            double p2f = cand2 / (double) sieveVectorBound;
            double p3f = cand3 / (double) sieveVectorBound;
            double p4f = cand4 / (double) sieveVectorBound;
            double expAll4 = p1f * p2f * p3f * p4f * sieveVectorBound;

            log();
            log("Block diag: thread", threadId,
                    " | baseX=", baseX,
                    " | smallHits=", smallHits,
                    " | bucketNodes=", bucketNodes);
            log("  cutoff≈", String.format("%.1f", cutoff * toBits),
                    " | candSUM=", candSUM,
                    " | candidatesForTD=", candidatesForTD,
                    " | gate=sum(r_i)≤cutoff");

            log("  Form-wise candidates (pre-TD):");
            log("    L1:", cand1, "  L2:", cand2, "  L3:", cand3, "  L4:", cand4);
            log("    ∩s:  L1∩L2:", cand12, "  L1∩L3:", cand13, "  L1∩L4:", cand14,
                    "  L2∩L3:", cand23, "  L2∩L4:", cand24, "  L3∩L4:", cand34);
            log("         L1∩L2∩L3:", cand123, "  L1∩L2∩L4:", cand124,
                    "  L1∩L3∩L4:", cand134, "  L2∩L3∩L4:", cand234,
                    "  ALL4:", cand1234, "  (E[ALL4]≈", String.format("%.3f", expAll4), " per block)");

            log("  Forms max bits in block [X∈{", MB.add(BigInteger.valueOf(xStart)), ", ",
                    MB.add(BigInteger.valueOf(xEnd)), "}]:");
            log("    L1=|sX-t|:", maxFormBits[0],
                    "  L2=|sX+t|:", maxFormBits[1],
                    "  L3=|fs²X-c|:", maxFormBits[2],
                    "  L4=|fs²X+c|:", maxFormBits[3],
                    "  | (sum max, intuition only):", (maxFormBits[0]+maxFormBits[1]+maxFormBits[2]+maxFormBits[3]));

            log("  Best-sum index details: idx=", bestI, " X=", Xbest);
            log("    L1:", bestFormBits[0], "  L2:", bestFormBits[1],
                    "  L3:", bestFormBits[2], "  L4:", bestFormBits[3],
                    "  | sum(bits):", bestSumBits,
                    "  | acc(sum bits):", String.format("%.1f", bestAccBits),
                    "  | rem(sum bits):", String.format("%.1f", bestRemBits));

            if (DEBUG_VERBOSE) {
                log("  Top", SHOW_TOP, "by product-style remainder (index->remBits):");
                for (int i = 0; i < SHOW_TOP; i++) {
                    if (bestIdx[i] < 0) break;
                    double rb = bestRem[i] * toBits;
                    log("   #", i + 1, ": idx=", bestIdx[i], " remBits≈", String.format("%.2f", rb));
                }
            }

            log("  TD (this block): B-smooth +", ws.debugSmooth,
                    ", 1-LP +", ws.debugLP1, ", 2-LP +", ws.debugLP2);
            ws.resetTDdebug();
            log();
        }
    }


    // ---- Expected ln*256 for each form at xOffset (block center) ----
    private int[] expectedFormLogsScaled(long xOffset) {
        BigInteger X = MB.add(BigInteger.valueOf(xOffset));
        BigInteger sX = s.multiply(X);
        BigInteger fs2X = f.multiply(s).multiply(s).multiply(X);
        double ln2 = Math.log(2.0);

        int b1 = sX.subtract(t).abs().bitLength();
        int b2 = sX.add(t).abs().bitLength();
        int b3 = fs2X.subtract(c).abs().bitLength();
        int b4 = fs2X.add(c).abs().bitLength();

        int e1 = (int) Math.round(b1 * ln2 * LOG_SCALE);
        int e2 = (int) Math.round(b2 * ln2 * LOG_SCALE);
        int e3 = (int) Math.round(b3 * ln2 * LOG_SCALE);
        int e4 = (int) Math.round(b4 * ln2 * LOG_SCALE);
        return new int[]{ e1, e2, e3, e4 };
    }

    private void trialDivideBucketed(long xOffset, int headNode, Work ws) {
        // X = M + x
        BigInteger X = MB.add(BigInteger.valueOf(xOffset));
        BigInteger sX   = s.multiply(X);
        BigInteger fs2X = f.multiply(s).multiply(s).multiply(X);

        // moving part only (product-level TD retained)
        BigInteger Qmov =
                sX.subtract(t).abs()
                        .multiply(sX.add(t).abs())
                        .multiply(fs2X.subtract(c).abs())
                        .multiply(fs2X.add(c).abs());

        if (Qmov.signum() == 0) return;

        BitSet bits = new BitSet(B_SMOOTH);
        BigInteger rem = Qmov;

        boolean foundFactor = false;

        // Step 1: small primes not bucketed
        for (int i = 0; i < primeBaseInt.length; i++) {
            int p = primeBaseInt[i];
            if (p >= SMALL_PRIME_LOG_ONLY_CUTOFF) break;

            if (p == 2) {
                int twoFactors = rem.getLowestSetBit();
                if (twoFactors > 0) {
                    rem = rem.shiftRight(twoFactors);
                    if ((twoFactors & 1) == 1) bits.set(i);
                    foundFactor = true;
                }
                continue;
            }

            int parity = 0;
            BigInteger pBI = primeBase[i];
            while (true) {
                BigInteger[] dr = rem.divideAndRemainder(pBI);
                if (dr[1].signum() == 0) {
                    rem = dr[0];
                    parity ^= 1;
                    foundFactor = true;
                } else break;
            }
            if (parity == 1) bits.set(i);
        }

        // Step 2: primes in bucket
        for (int node = headNode; node != -1; node = ws.next[node]) {
            int i = ws.who[node];
            int parity = 0;
            BigInteger pBI = primeBase[i];
            while (true) {
                BigInteger[] dr = rem.divideAndRemainder(pBI);
                if (dr[1].signum() == 0) {
                    rem = dr[0];
                    parity ^= 1;
                    foundFactor = true;
                } else break;
            }
            if (parity == 1) bits.set(i);
        }

        // B-smooth?
        if (rem.equals(BigInteger.ONE)) {
            ws.localBSmooth.add(new VectorData(bits, xOffset));
            ws.localBSmoothCount++;
            ws.debugSmooth++;
            return;
        }

        // No division? sample sometimes
        if (!foundFactor && DEBUG_VERBOSE && (rnd.nextInt(2048) == 0)) {
            log("  TD sample: rem(no-division) bits≈", rem.bitLength(), " | xOffset=", xOffset);
        }

        // LP logic
        if (rem.compareTo(largePrimeBoundBI) > 0) return;

        if (rem.isProbablePrime(20)) {
            ws.localBigPrimes.add(new Pair(rem.longValue(), new VectorData(bits, xOffset)));
            ws.debugLP1++;
            return;
        }

        long r = rem.longValue();
        PairLong pq = factorSemiprimeLE1e10(r);
        if (pq == null) return;
        if (!BigInteger.valueOf(pq.a).isProbablePrime(20)) return;
        if (!BigInteger.valueOf(pq.b).isProbablePrime(20)) return;

        VectorData cur = new VectorData(bits, xOffset);
        VectorData merged = bigPrimePairs.add(pq.a, pq.b, cur);
        if (merged != null) {
            ws.localBSmooth.add(merged);
            ws.localBSmoothCount++;
            ws.debugLP2++;
        }
    }

    private PairLong factorSemiprimeLE1e10(long n) {
        if ((n & 1L) == 0) return new PairLong(2L, n >>> 1);
        if (n % 3L == 0) return new PairLong(3L, n / 3L);
        if (n % 5L == 0) return new PairLong(5L, n / 5L);
        for (long d = 7; d * d <= n && d <= 100_000L; d += 2) {
            if (n % d == 0) return new PairLong(d, n / d);
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
            if (solution.isEmpty()) continue;
            log("Testing solution", (i + 1) + "/" + solutions.size());
            if (testSolution(solution)) return true;
        }
        log("no luck");
        return false;
    }

    private boolean isReadyToBeSolved() {
        return bSmoothVectors.size() + bigPrimesList.getPrimesFound() >= B_SMOOTH + 10;
    }

    // A(X) = f^2 * s^2 * (s X - t)(s X + t)  (mod N) — modular sqrt of Q5(X)
    private BigInteger A_of(BigInteger X) {
        BigInteger sX = s.multiply(X);
        BigInteger term = sX.subtract(t).multiply(sX.add(t));
        return f.multiply(f).multiply(s).multiply(s).multiply(term).mod(N);
    }

    // Q5(X) = (f^2 s^2) * moving part (integer, for Y-side)
    private BigInteger Q5_of(BigInteger X) {
        BigInteger sX   = s.multiply(X);
        BigInteger fs2X = f.multiply(s).multiply(s).multiply(X);
        BigInteger mov  = sX.subtract(t).multiply(sX.add(t))
                .multiply(fs2X.subtract(c))
                .multiply(fs2X.add(c));
        return f.multiply(f).multiply(s).multiply(s).multiply(mov);
    }

    private boolean testSolution(ArrayList<VectorData> solutionVector) {
        BigInteger one = BigInteger.ONE;
        BigInteger Xside = one;
        BigInteger Yprod = one;

        for (VectorData vd : solutionVector) {
            BigInteger ax, qx;
            if (vd.x != null) {
                ax = vd.x;     // already A-product mod N
                qx = vd.y;     // already Q5-product (integer)
            } else {
                BigInteger X = MB.add(BigInteger.valueOf(vd.position));
                ax = A_of(X);
                qx = Q5_of(X);
            }
            Xside = Xside.multiply(ax).mod(N);
            Yprod = Yprod.multiply(qx);
        }

        BigInteger Yside = SqrRoot.bigIntSqRootFloor(Yprod.abs());

        BigInteger g = N.gcd(Xside.subtract(Yside).abs());
        if (!g.equals(one) && !g.equals(N)) {
            log("Solved");
            log(g);
            return true;
        }
        g = N.gcd(Xside.add(Yside).abs());
        if (!g.equals(one) && !g.equals(N)) {
            log("Solved");
            log(g);
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

    // ===== Debug helpers =====

    private void calibrateAndLogOnce() {
        // Center X = M
        BigInteger X = MB;
        BigInteger sX   = s.multiply(X);
        BigInteger fs2X = f.multiply(s).multiply(s).multiply(X);

        int l1 = sX.subtract(t).abs().bitLength();
        int l2 = sX.add(t).abs().bitLength();
        int l3 = fs2X.subtract(c).abs().bitLength();
        int l4 = fs2X.add(c).abs().bitLength();

        double targetBits = Math.max(1, N.bitLength() / 5.0);

        log("Calibration (at center X=M):");
        log("  L1(sX-t)=", l1, " L2(sX+t)=", l2, " L3(fs^2X-c)=", l3, " L4(fs^2X+c)=", l4);
        log("  moving-sum bits≈", l1 + l2 + l3 + l4, " (target≈", String.format("%.1f", targetBits * 4), " for 4 linears)");
        double cutoffBits = (BIG_PRIME_CUTOFF_SCALED + LOG_EPS) / (Math.log(2.0) * LOG_SCALE);
        log("  large-prime cutoff (bits)≈", String.format("%.2f", cutoffBits), " (this is log2(B_max^2))");

        for (int i = 0; i < 3; i++) {
            long dx = (rnd.nextInt(2) == 0 ? -1 : 1) * (1L << (i * 4 + 6)); // ±64, ±1k, ±16k
            BigInteger Xs = MB.add(BigInteger.valueOf(dx));
            int msum = movingSumBits(Xs);
            log("  sample X=M+", dx, ": moving-sum bits≈", msum);
        }
    }

    private void logCenteringFeasibility() {
        int bt = t.abs().bitLength();
        int bc = c.abs().bitLength();
        int bs = s.abs().bitLength();
        int bf = f.abs().bitLength();

        int rhs = bc - bt;
        int lhs = bf + bs;

        log("Centering diagnostic (bits):");
        log("  bits(t)=", bt, " bits(c)=", bc, " bits(s)=", bs, " bits(f)=", bf);
        log("  Need  : log2(f)+log2(s) ≈ bits(c)-bits(t) =", rhs);
        log("  Have  : log2(f)+log2(s) =", lhs, "  (mismatch=", (lhs - rhs), ")");
        if (rhs < 0) {
            log("  => IMPOSSIBLE TO CENTER: c is ~", -rhs, " bits smaller than t; one pair stays large.");
        } else if (lhs < rhs) {
            log("  => Under-sloped: need ~", (rhs - lhs), " more bits in f*s to align both pairs.");
        } else {
            log("  => Feasible in principle; fine centering depends on exact values and M.");
        }
    }

    private int movingSumBits(BigInteger X) {
        BigInteger sX   = s.multiply(X);
        BigInteger fs2X = f.multiply(s).multiply(s).multiply(X);
        return sX.subtract(t).abs().bitLength()
                + sX.add(t).abs().bitLength()
                + fs2X.subtract(c).abs().bitLength()
                + fs2X.add(c).abs().bitLength();
    }

    private void wheelSanity(SingleWheel[][] wheels) {
        int[] testPs = {3,5,7,11};
        long baseX = 0;
        for (int p : testPs) {
            if (p >= primeBaseInt[B_SMOOTH - 1]) break;
            int idx = -1;
            for (int i = 0; i < primeBaseInt.length; i++) if (primeBaseInt[i] == p) { idx = i; break; }
            if (idx < 0) continue;
            SingleWheel w = wheels[idx][0]; // check L1
            if (w == null) { log("Wheel sanity: p=", p, " (L1) has no residues (OK if never divides)"); continue; }
            w.rebase(baseX, sieveVectorBound);
            int count = 0;
            while (w.hasNext()) { w.next(); count++; }
            log("Wheel sanity: p=", p, " hits in first block (L1):", count);
        }
    }

    // ===== Work buffers =====
    private static final class Work {
        final int[] accL1, accL2, accL3, accL4; // per-form scaled log sums
        final int[] head;                       // product-level bucket head
        final ArrayList<VectorData> localBSmooth = new ArrayList<>(256);
        final ArrayList<Pair>       localBigPrimes = new ArrayList<>(256);
        int[] who;                    // prime index for node
        int[] next;                   // next node pointer
        int ptr = 0;
        int localBSmoothCount = 0;

        // TD debug counters
        int debugSmooth = 0;
        int debugLP1 = 0;
        int debugLP2 = 0;

        Work(int n, int bucketCap) {
            this.accL1 = new int[n];
            this.accL2 = new int[n];
            this.accL3 = new int[n];
            this.accL4 = new int[n];
            this.head = new int[n];
            this.who = new int[bucketCap];
            this.next = new int[bucketCap];
        }

        void growBuckets() {
            int newCap = who.length * BUCKET_GROWTH;
            who  = Arrays.copyOf(who,  newCap);
            next = Arrays.copyOf(next, newCap);
        }

        void resetTDdebug() { debugSmooth = debugLP1 = debugLP2 = 0; }
    }

    private record PairLong(long a, long b) {}
    private record Pair(long prime, VectorData data) {}

    // ===== Single linear wheel: solutions of a*(M+x) ≡ b (mod p) in x =====
    private static final class SingleWheel {
        final int p;
        final int residue;   // unique residue for x modulo p, or -1 if no solution, or -2 for "all x"
        int baseMod;
        int bound;

        int cur;
        int left;

        private SingleWheel(int p, int residue) {
            this.p = p; this.residue = residue;
        }

        static SingleWheel build(int p, BigInteger a, BigInteger b, BigInteger MB) {
            if (p <= 1) return null;
            int aMod = mod(a, p);
            int bMod = mod(b, p);
            int mMod = mod(MB, p);

            if (aMod == 0) {
                // equation: 0 * (M+x) ≡ b (mod p)  => b ≡ 0  => either all x or none
                return new SingleWheel(p, (bMod == 0) ? -2 : -1);
            }
            int invA = invMod(aMod, p);
            if (invA == -1) return new SingleWheel(p, -1); // shouldn't happen for prime p

            // a*(M+x) ≡ b  => a*x ≡ b - a*M
            int rhs = bMod - (int)((long)aMod * mMod % p);
            rhs %= p; if (rhs < 0) rhs += p;
            int r = (int)((long)rhs * invA % p);
            return new SingleWheel(p, r);
        }

        void rebase(long baseX, int blockSize) {
            this.baseMod = mod(baseX, p);
            this.bound = blockSize;

            if (residue == -1) { left = 0; cur = 0; return; }
            if (residue == -2) {
                // all x are solutions
                cur = 0;
                left = blockSize;
                return;
            }
            int idx = residue - baseMod; if (idx < 0) idx += p;
            cur = idx;

            int count = 0;
            if (idx < blockSize) count = 1 + (blockSize - 1 - idx) / p;
            left = count;
        }

        boolean hasNext() { return left > 0; }
        int next() {
            if (residue == -2) { // all x in block
                int out = cur;
                cur++;
                left--;
                return out;
            }
            int out = cur;
            cur += p;
            left--;
            return out;
        }

        private static int mod(BigInteger x, int p) {
            int m = x.mod(BigInteger.valueOf(p)).intValue();
            if (m < 0) m += p;
            return m;
        }
        private static int mod(long x, int p) { int m = (int)(x % p); if (m < 0) m += p; return m; }

        private static int invMod(int a, int p) {
            int t = 0, newT = 1;
            int r = p, newR = a;
            while (newR != 0) {
                int q = r / newR;
                int tmpT = t - q * newT; t = newT; newT = tmpT;
                int tmpR = r - q * newR; r = newR; newR = tmpR;
            }
            if (r != 1) return -1;
            if (t < 0) t += p;
            return t;
        }
    }

    // Bits of each linear at a specific X (moving part only)
    private int[] formBitsAtX(BigInteger X) {
        BigInteger sX   = s.multiply(X);
        BigInteger fs2X = f.multiply(s).multiply(s).multiply(X);
        int b1 = sX.subtract(t).abs().bitLength();   // |sX - t|
        int b2 = sX.add(t).abs().bitLength();        // |sX + t|
        int b3 = fs2X.subtract(c).abs().bitLength(); // |fs^2 X - c|
        int b4 = fs2X.add(c).abs().bitLength();      // |fs^2 X + c|
        return new int[]{ b1, b2, b3, b4 };
    }

    // Max bits per form over the whole block X-range (achieved at an endpoint for linear |aX±b|)
    private int[] formMaxBitsForBlock(long xStart, long xEnd) {
        BigInteger X0 = MB.add(BigInteger.valueOf(xStart));
        BigInteger X1 = MB.add(BigInteger.valueOf(xEnd));
        int[] a = formBitsAtX(X0);
        int[] b = formBitsAtX(X1);
        return new int[]{
                Math.max(a[0], b[0]),
                Math.max(a[1], b[1]),
                Math.max(a[2], b[2]),
                Math.max(a[3], b[3])
        };
    }
}
