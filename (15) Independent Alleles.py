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
