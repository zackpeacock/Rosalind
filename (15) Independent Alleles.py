# Given: Two positive integers k (k ≤ 7) and N (N ≤ 2^k). In this problem, we begin with Tom, who in the 0th generation has genotype Aa Bb. 
# Tom has two children in the 1st generation, each of whom has two children, and so on. Each organism always mates with an organism having genotype Aa Bb.

# Return: The probability that at least N Aa Bb organisms will belong to the k-th generation of Tom's family tree (don't count the Aa Bb mates at each level). 
# Assume that Mendel's second law holds for the factors.


# Assign given values for k and N to variables
import math
file = open('rosalind_lia.txt', 'r')
givenIntegers = file.read().strip().split(' ')

k = int(givenIntegers[0])
N = int(givenIntegers[1])
P = 2 ** k

# Probability of heterozygous offspring if at least 1 parent is heterozygous = 0.25
# Probability of N heterozygotes in total population P = [P! / N!(P - N)!] * (0.25 ^ i) * [0.75 ^ (P - i)]
# Sum probabilities of producing N, N + 1, N + 2 ... P heterozygotes
probability = 0
for i in range(N, P + 1):
    single_gen_prob = (math.factorial(P) / (math.factorial(i) * math.factorial(P - i))) * (0.25 ** i) * (0.75 ** (P - i))
    probability += single_gen_prob

print(round(probability, 3))
