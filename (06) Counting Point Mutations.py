# Given: Two DNA strings of equal length (not exceeding 1 kbp).

# Return: The Hamming distance between the two strings.


file = open('rosalind_hamm.txt', 'r')

# Assign sequences to be compared to variables and initialize mutation counter
rawSequences = file.read().split('\n')
originalSequence = rawSequences[0]
mutatedSequence = rawSequences[1]
mutations = 0

# Compute number of mutations observed between sequences
for i, value in enumerate(originalSequence):
    if originalSequence[i] == mutatedSequence[i]:
        pass
    else:
        mutations += 1

print(mutations)
