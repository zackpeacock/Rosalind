# Given: A positive integer n ≤ 7.

# Return: The total number of permutations of length n, followed by a list of all such permutations (in any order).


from itertools import permutations
import math

# Assign given integer to variable n
file = open('rosalind_perm.txt', 'r')
given_integer = file.read().split()
n = int(given_integer[0])

# Create list of permutations of n
perms = list(permutations(range(1, n + 1)))

# Format list of permutations
formatted_perms = []
for lst in perms:
    sublist = []
    for i in lst:
        sublist.append(str(i))
    formatted_perms.append(sublist)

# Print number of permutations and all possible permutations in proper format
print(math.factorial(n))
for lst in formatted_perms:
    print(' '.join(lst))
