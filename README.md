# Java Quadratic Sieve

This is Quadratic sieve implementation in Java, you can find a lot of info about the algorithm at [wiki](https://en.wikipedia.org/wiki/Quadratic_sieve), so I just assume
you are familiar with the basics of the algorithm and jump right in to the implementation.
 
# Implementation 

 - First I start by building the prime base. I am using a predefined B_BOUNDS value to determinate how many primes I want to include in it. I iterate over all the numbers starting from 2, and using the [primality test](https://en.wikipedia.org/wiki/Primality_test) to check if they are prime(BigInteger.nextProbablePrime). I am using [Legendre symbol](https://en.wikipedia.org/wiki/Legendre_symbol) to include only those primes that are in quadratic residue mod N.
 - Next I am building wheels, each wheel holds one prime, it tells what is the next number that the prime will divide. I use them during the sieving loop. In order to find the first tow values of each prime I use [Tonelli–Shanks](https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm) algorithm, he called it ressol. Also the wheel calculates the prime log. This is in order to avoid unnecessary division during the sieve loop.
 - Now I star the sieve loop, I start checking values from the ceil root of N(the values that I am trying to factor) and define **sieveVectorBound** as the size of the biggest prime in the factory base, this is how many numbers I will process over each loop. I choose this size so each wheel will be used attlist once.
   - First I calculate the base log, it's the size ceil root of N plus current **position**, position stats from 0 and increased by sieveVectorBound each loop.
   - I create 2 arrays of doubles, **logs** and **trueLogs** with the size of sieveVectorBound. 
   - I create another array of **VectorData** with the size of sieveVectorBound. VectorData contains information about the position where it been found, and the b-smooth vector already transformed to boolean, I use **BitSet** for it. Position is used to calculate **x** and **y** such that **x^2 - N = y**, it's important for the last step of the algorithm.
   - Now I iterate over all the wheels and sum all their logs in the logs array. I do it with the help of tow important methods that wheel got:
     - testMove(long *limit*) it checks if the wheel already reached a *limit* position
     - nextLog() it increase the wheel current position and return the log for that position(At this point its the same log all the time, but in early implementations I been returning the log multiplied by the biggest power of that prime that could divide that position)   
   - Next I return the wheels to the position they been at when this loop started and make another iteration over them all
   - For each wheel I check if it's corresponding log[index] is bigger than the log of biggest prime in the prime base power tow. If it's bigger I skip to the next wheel but if it's not, this number is aether a prime or divided by one of the prime base primes.
   - So if not, I find the actual log of that position by taking it from the trueLogs array or calculating and saving it there if it's the first time. I compare the actual log to the log I calculated if its equal withing a small error of **0.0000001**, I know that this position is fully factored by the prime base and I add its VectorData to bSmoothVectors VectorData array to be used latter.
   - if its not equal than it's a prime, so I add its VectorData to **bigPrimesList**.   
   - We are still in the if case of not ignoring the wheel, we added the VectorData to the right list and now I turn on the bit of VectorData b-smooth vector at the current index of the wheels array
 - When enough big primes and b-smooth vectors been found, I construct and solve the matrix using [gaussian elimination](https://en.wikipedia.org/wiki/Gaussian_elimination), I know it can be improved by using algorithms to find the null space, but I am using a very efficient implementation with bit calculation, so it's not the bottle neck of the algorithm
 - I make small optimization in the bigPrimeList, I only take those big primes that been found more than once, and if the big prime been found just twice, I xor, those two vectors to one, as it's the only way it will be part of the solution.
      
 After solving the matrix, I extract the solution and it's pretty much it. As of version 2.1 I can factor 140 bits numbers(40 digits) in half minute. The matrix solver can handle matrices of the size 100,000 in 15 minutes.
 
# How can you help?
 
This code needs better documentation, the implementation description can be improved as well, I did it because I love math.
I will appreciate any bug reports, comments and new ideas.    

