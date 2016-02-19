# quadratic-sieve
Quadratic sieve implementation in Java

<B>As of version 2.1 it can factor 40 digits number in 30 sec and solve matrixed with the size of 100,000x100,000 in 15 min.</B>

This is java implemintation of the [Quadratic Sieve](https://en.wikipedia.org/wiki/Quadratic_sieve) algorithm.

The major improvement over the regular sieve is the prime wheels, those allows me to avoid the computision of values that are not part of B-Smooth vector.
Perhaps I should start by showing how this algorithm is actually working.

 - First the algorithm find a predefined number of primes, starting from 1 up to some const value B_SMOOTH, during the search it test to make sure that those primes are in rootIn quadratic residues of X^2 - N
 - Than it create a predefined number of VectorNumber objects and storing them in vectorNumbers array, the size of the vector represented by the const SIEVE_VECTOR_BOUND
 - A wheel is created for each prime in the prime base. The wheel can tell what would be the next value that its prime(Each wheel got one prime that it's hosting) will divide with mod 0.
 - Now the search for the B-Smooth numbers is started, and each step looks as fallow
   - First the vector version is increased, this reset all the values in the vectores array in simple O(1) operation
   - The current position is increased by SIEVE_VECTOR_BOUND
   - Each wheel is rolling, until it reached the new limit, and during its roll it divide all the values in the vectorNumbers, so if one of the values become 1 a B-Smooth is extracted from it.
 - After there are more than B_SMOOTH values been found. An attempt to solve the matrix is performed
 - The algorithm terminates when non trivial factor is found
 
 So far I been able to factor values up to 30 digits length.
 
 This algorithm can be improved by adding b-smooth values that got one [large prime](https://en.wikipedia.org/wiki/Quadratic_sieve#One_large_prime). Another point of improvement is its memory consumption. 
 It keeps the vectors of each value of SIEVE_VECTOR_BOUND in memory, so while it make it a bit faster in the last step, it slows it during the search process and limiting the prime base size.
 
 So while there are still a lot to improve this algorithm simply works. And it's doing a great job! 
 I hope you enjoy it just the same way I do.
__________________________________________

V1.01
----------

Improved the memory performance by moving the in search vector building part, where I break b-smooth values to vector. I moved it to after the search. So it slow down the performance since it's extra step on the way, but it allows to factor bigger numbers since memory is no longer a limit.

V1.5
---------

 - Fixed bug with wheels, Sometimes wheel 2 been faster than wheel 1 and it caused to miss some factors.
 - Added Stefan Buettcher implementation of Tonelli–Shanks algorithm to speed up the wheel initialization. 
   I made some small optimizations to his code, but there are parts that are pretty messy, like MathUtils.v_ method that I don't fully understand yet.

At this point the bottleneck is in solving the matrix, it turns out that the current implementation doing it in O(N^4) so it's taking hours to solve 10K vectors.

V1.6
----------

This version is all about matrixes, I added two new versions of it. The <b>HashMatrix</b>, it supposed to take less memory and have better performance as it doesn't iterate over the zeroes in the matrix, but in practice the hashing operation is very expensive and the new BitMatrix perform much better, both in terms of memory and performance. With <b>BitMatrix</b> I been able to solve 5K B bound with about 5K rows in minutes. 

Now there is a new challenge, after the solution been calculated there is need to multiply the solution rows and then square them in order to get the Y value for gcd, for X is a bit easier as I keep it in original form, so no need to square it. 
For values over 150 bit length, it takes more then minutes to squere. 

V2.0
----------

In this version I discovered the power of logarithms. Instead of dividing you can sum the logs in order to tell if some value is divisible by others.

I also found the bottleneck of bitMatrix, it turns out that 90% of the calculation time went on initializing, so I completely throw it away. Instead BigInteger I moved to BitSet, it supposed to be much faster how ever I don’t use the BitSet.nextSet method as it will prevent me from using xor, so it makes them pretty much the same in performance, how ever the BitSet allows me to initialize it with maximum bit size, and this makes a big difference, it will used less memory.

So now I build the vectors during the searching loop process and when I come to calculation part, I already got the bitMatrix almost initialized. This allows me to solve matrixes with the size of 100,000 rows in about 15 min.

The other problem that I had is the solution extraction, I used babilon method for root squaring and it took hours. Now I did a very small improvement to it, I changed the first guess of the square to be 2, power half the bit size of the number. This reduced the calculation time from hours to less than a millisecond, I been amazed to see this, how ever I suspect that the difference will not be as great with numbers that don’t got an integer square. But it’s not the case here, and the extraction works now in seconds. It yet might be a problem for numbers bigger than 60 decimal digits.

Now the bottleneck is finally where it should be, in the process of finding the B-Smoothe vectors.

The next tasks are:

 - Write a better description for the algorithm, as at this point I probably the only one who can understand it’s description…
 - Add Big primes optimization

V2.1
----------
Added big primes optimization. It turns out that the wheels did a lot of work while predicting the exact power that the prime can be divided by, so I removed this future and it become faster. 



