# Given: A collection of k (k ≤ 100) DNA strings of length at most 1 kbp each in FASTA format.

# Return: A longest common substring of the collection.


# Create dictionary of sequence labels with sequences
file = open('rosalind_lcsm.txt', 'r')
sequencesDict = {}
sequenceLabel = ''
for line in file.readlines():
    if '>' in line:
        sequenceLabel = line.rstrip()
        sequencesDict[sequenceLabel] = ''
    else:
        sequencesDict[sequenceLabel] += line.rstrip()

# Create list of sequences
sequences = []
for value in sequencesDict.values():
    sequences.append(value)

# Sort sequences from shortest to longest, identify shortest as sequence to compare to all others
sorted_sequences = sorted(sequences, key=len)
shortest_sequence = sorted_sequences[0]
other_sequences = sorted_sequences[1:]

# Create list of all possible substrings of shortest sequence
possible_substrings = [shortest_sequence[i:j] for i in range(len(shortest_sequence))
                       for j in range(i + 1, len(shortest_sequence) + 1)]

# Remove duplicate substrings and sort list from longest to shortest
unique_substrings = list(set(possible_substrings))
sorted_substrings = sorted(unique_substrings, key=len, reverse=True)

# Identify longest substring (motif) present in all other sequences
motif = ''
for substring in sorted_substrings:
    for sequence in sequences:
        if substring in sequence:
            match = True
        else:
            match = False
            break
    if match and len(substring) > len(motif):
        motif = substring

print(motif)
