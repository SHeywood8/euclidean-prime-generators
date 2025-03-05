import math
import numpy as np

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
    '''Simple prime factorisation function through brute force division checks'''
    factors = []
    for p in eratosthanes(n):
        while n % p == 0:
            factors.append(p)
            n //= p
    return factors

def divisors(n):
    '''Simple divisor function through brute force division checks'''
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