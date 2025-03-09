import math
import numpy as np
import random

def eratosthanes(n):
    '''Simple prime sieve function for calculating all primes up to n'''
    n = int(n)
    sieve = [i for i in range(0,n+1)]
    sieve[1] = 0
    for x in sieve:
        if x != 0:
            for i in range(2,n+1):
                if x*i > n:
                    break
                sieve[x*i] = 0
    return [x for x in sieve if x != 0]

def factor(n):
    '''Simple prime factorisation function through brute force division checks
    Only intended for use on small n'''
    factors = []
    for p in eratosthanes(n):
        while n % p == 0:
            factors.append(p)
            n //= p
    return factors

def divisors(n):
    '''Simple divisor function through brute force division checks
    Only intended for use on small n'''
    divs = []
    for d in range(1,n):
        if n % d == 0:
            divs.append(d)
    return divs

def Jacobi(a,n): # Function from https://literateprograms.org/jacobi_symbol__python_.html
    '''Returns the calculated Jacobi symbol (a|n) 
    Uses properties of the symbol and Quadratic Reciprocity'''
    if a == 0:
        return 0
    if a == 1:
        return 1
    if a == 2:
        n8 = n % 8
        if n8 in [3,5]:
            return -1
        else:
            return 1
    if a % 2 == 0:
        return Jacobi(2,n)*Jacobi(a//2,n)
    if a >= n:
        return Jacobi(a % n, n)
    if a % 4 == 3 and n % 4 == 3:
        return (-1)*Jacobi(n,a)
    else:
        return Jacobi(n,a)

def fermat(a : int, d : int, N : int): # Book page 88
    '''Calculates a^d mod N'''
    prod = 1
    a2j = a
    while d>0:
        if d % 2 == 1:
            prod *= a2j
            prod %= N
        d //= 2
        a2j **= 2
        a2j %= N
    return prod

def fermat_test(N : int, iter = 10, baselim = 313):
    '''Performs a fermat compositeness test with multiple prime bases
     to try to rule out if N is prime. Output is whether or not N is
     (probably) prime.'''
    primes = eratosthanes(baselim)
    for p in primes:
        if N % p == 0 and N != p:
            return False
    for i in range(iter):
        if math.gcd(primes[i],N) == 1:
            if fermat(primes[i],N-1,N) != 1:
                return False
    else:
        return True

def modular_pow(base, exponent, modulus): # Function from https://www.geeksforgeeks.org/pollards-rho-algorithm-prime-factorization/
    '''Efficiently computes large powers base**exponent mod modulus'''
    result = 1  # initialize result 
    while (exponent > 0):
        if (exponent & 1): # if y is odd, multiply base with result 
            result = (result * base) % modulus
        exponent = exponent >> 1 # exponent = exponent/2 
        base = (base * base) % modulus # base = base * base 
    return result

def PollardRho(n): # Function from https://www.geeksforgeeks.org/pollards-rho-algorithm-prime-factorization/
    '''Implementation of the Pollard Rho algorithm to find factors'''
    if (n == 1): # no prime divisor for 1 
        return n
    if (n % 2 == 0): # even number means one of the divisors is 2 
        return 2
    x = (random.randint(0, 2) % (n - 2)) # we will pick from the range [2, N) 
    y = x
    # the constant in f(x).
    # Algorithm can be re-run with a different c
    # if it throws failure for a composite. 
    c = (random.randint(0, 1) % (n - 1))
    # Initialize candidate divisor (or result) 
    d = 1
    # until the prime factor isn't obtained.
    # If n is prime, return n 
    while (d == 1):
        # Tortoise Move: x(i+1) = f(x(i)) 
        x = (modular_pow(x, 2, n) + c + n)%n
        # Hare Move: y(i+1) = f(f(y(i))) 
        y = (modular_pow(y, 2, n) + c + n)%n
        y = (modular_pow(y, 2, n) + c + n)%n
        # check gcd of |x-y| and n 
        d = math.gcd(abs(x - y), n)
        # retry if the algorithm fails to find prime factor
        # with chosen x and c 
        if (d == n):
            return PollardRho(n)
    return d

def factor2(n):
    '''Function that finds the prime factorisation of n.
    First brute force checks all primes less than 10000
    Then checks if remaining number is composite
    Then repeatedly uses PollardRho algorithm to find larger factors'''
    factors = []
    for p in eratosthanes(10000):
        while n % p == 0:
            n //= p
            factors.append(p)
    while n > 1:
        if not fermat_test(n):
            a = PollardRho(n)
            factors.append(a)
            n //= a
        else:
            factors.append(n)
            return factors
    return factors