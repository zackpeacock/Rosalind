# Given: Two DNA strings s and t (each of length at most 1 kbp).

# Return: All locations of t as a substring of s.


file = open('rosalind_subs.txt', 'r')

# Assign sequence and subsequence to variables
rawSequences = file.read().split('\n')
sequence = rawSequences[0]
subsequence = rawSequences[1]

# Create list of locations of subsequence
locations = []
n = len(subsequence)
for i, value in enumerate(sequence):
    if sequence[i:i+n] == subsequence:
        locations.append(str(i+1))

print(' '.join(locations))
