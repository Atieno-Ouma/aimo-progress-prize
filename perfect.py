import math

# Task 1: Check if a number is a perfect number
def is_perfect_number(n):
    divisors_sum = sum(i for i in range(1, n) if n % i == 0)
    return divisors_sum == n

# Task 2: Generate Pascal's Triangle up to n rows
def generate_pascals_triangle(n):
    triangle = [[1]]
    for i in range(1, n):
        row = [1]
        for j in range(1, i):
            row.append(triangle[i-1][j-1] + triangle[i-1][j])
        row.append(1)
        triangle.append(row)
    return triangle

# Task 3: Check if a number is a Carmichael number (A type of pseudoprime)
def is_carmichael_number(n):
    if n < 2:
        return False
    for a in range(2, n):
        if math.gcd(a, n) == 1 and pow(a, n-1, n) != 1:
            return False
    return True

# Task 4: Find the nth triangular number
def triangular_number(n):
    return n * (n + 1) // 2

# Task 5: Find the least common multiple (LCM) of two numbers
def lcm(a, b):
    return abs(a * b) // math.gcd(a, b)

# Task 6: Solve a quadratic equation (ax^2 + bx + c = 0)
def solve_quadratic(a, b, c):
    discriminant = b**2 - 4*a*c
    if discriminant > 0:
        root1 = (-b + math.sqrt(discriminant)) / (2*a)
        root2 = (-b - math.sqrt(discriminant)) / (2*a)
        return root1, root2
    elif discriminant == 0:
        root = -b / (2*a)
        return root,
    else:
        return "No real roots"

# Task 7: Generate prime numbers up to a given limit using the Sieve of Eratosthenes
def sieve_of_eratosthenes(limit):
    primes = [True] * (limit + 1)
    primes[0], primes[1] = False, False
    for i in range(2, int(math.sqrt(limit)) + 1):
        if primes[i]:
            for j in range(i * i, limit + 1, i):
                primes[j] = False
    return [i for i in range(limit + 1) if primes[i]]

# Task 8: Compute the binomial coefficient (n choose k)
def binomial_coefficient(n, k):
    if k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# Task 9: Solve a system of linear equations using matrix inversion (3x3 matrix)
import numpy as np

def solve_system_with_matrix_inversion(a, b):
    try:
        a_inv = np.linalg.inv(a)
        solution = np.dot(a_inv, b)
        return solution
    except np.linalg.LinAlgError:
        return "The matrix is singular, no solution."

# Task 10: Compute the power of a number (x^n)
def power(x, n):

