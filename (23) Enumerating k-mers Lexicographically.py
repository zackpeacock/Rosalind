# Given: A collection of at most 10 symbols defining an ordered alphabet, and a positive integer n (n ≤ 10).

# Return: All strings of length n that can be formed from the alphabet, ordered lexicographically (use the standard order of symbols in the English alphabet).


import itertools

# Create and format alphabet list from .txt file, assign string length to variable
file = open('rosalind_lexf.txt', 'r')
alphabet = file.read().strip().split(' ')
length = alphabet[-1][-1]
alphabet[-1] = alphabet[-1][0]

# Create list of all possible permutations of given string length, including those in which letters repeat
# Python automatically sorts lexicographically
permutations = []
for i in itertools.product(alphabet, repeat=int(length)):
    permutations.append(i)

# Convert permutations to string form
strings = []
for i in permutations:
    string = ''.join(i)
    strings.append(string)

# Print individual strings
for string in strings:
    print(string)
