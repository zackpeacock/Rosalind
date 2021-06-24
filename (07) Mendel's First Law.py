# Given: Three positive integers k, m, and n, representing a population containing k + m + n organisms: k individuals are homozygous dominant for a factor, 
# m are heterozygous, and n are homozygous recessive.

# Return: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele 
# (and thus displaying the dominant phenotype). Assume that any two organisms can mate.


integers = open('rosalind_iprb.txt', 'r')
for string in integers:
    integerList = string.strip().split(' ')

# Assign numbers of homozygous dominant (k), heterozygous (m), and homozygous recessive (n) individuals to variables
k = int(integerList[0])
m = int(integerList[1])
n = int(integerList[2])

# t = population total
t = k + m + n

# Probabilities k, m, or n is first selection
probk1 = k / t
probm1 = m / t
probn1 = n / t

# Probabilities k, m, or n is second selection if not first selection
probk2 = k / (t - 1)
probm2 = m / (t - 1)
probn2 = n / (t - 1)

# Probabilities k, m, or n is second selection if also first selection
probk3 = (k - 1) / (t - 1)
probm3 = (m - 1) / (t - 1)
probn3 = (n - 1) / (t - 1)

# Probability of two recessive alleles
probRec = (0.25 * (probm1 * probm3)) + (0.5 * (probm1 * probn2)) + (probn1 * probn3) + (0.5 * (probn1 * probm2))

# Probability of one dominant allele
probDom = 1 - probRec

print(round(probDom, 5))
