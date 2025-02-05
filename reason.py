import math

# Task 1: Check if a matrix is symmetric
def is_symmetric(matrix):
    n = len(matrix)
    for i in range(n):
        for j in range(i, n):  # Only need to check the upper triangle
            if matrix[i][j] != matrix[j][i]:
                return False
    return True

# Task 2: Check if a number is prime
def is_prime(n):
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

# Task 3: Fibonacci sequence up to n-th term
def fibonacci(n):
    fib_seq = [0, 1]
    while len(fib_seq) <= n:
        fib_seq.append(fib_seq[-1] + fib_seq[-2])
    return fib_seq[:n]

# Task 4: Check if a number is a perfect square
def is_perfect_square(n):
    return n == int(math.sqrt(n)) ** 2

# Task 5: Calculate the greatest common divisor (GCD) of two numbers
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# Task 6: Solve a system of linear equations using Gaussian elimination (2x2 matrix for simplicity)
def solve_linear_system(a, b):
    # Solve ax + by = c and dx + ey = f for x and y
    # a, b = coefficients, c, d = constant terms
    det = a[0][0]*a[1][1] - a[0][1]*a[1][0]

