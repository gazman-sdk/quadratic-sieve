# quadratic-sieve
Quadratic sieve implementation in Java

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

 - Fixed bug with wheels, that sometimes wheel 2 been faster than wheel 1 and it caused to miss some factors.
 - Added Stefan Buettcher implementation of Tonelli–Shanks algorithm to speed up the wheel initialization. 
   I made some small optimizations to his code, but there are parts that are pretty messy, like MathUtils.v_ method that I don't fully understand yet.

At this point the bottleneck is in solving the matrix, it turns out that the current implementation doing it in O(N^4) so it's taking hours to solve 10K vectors.





