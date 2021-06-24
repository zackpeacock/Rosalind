# Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.

# Return: A consensus string and profile matrix for the collection.


import numpy as np
import pandas

file = open('rosalind_cons.txt', 'r')

# Create dictionary of sequence labels with sequences
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

# Create keys for nested dictionary (nucleotide_frequency) from length of sequence
length = []
def sequenceLength(first_nucleotide_index, last_nucleotide_index):
    while first_nucleotide_index < last_nucleotide_index:
        length.append(first_nucleotide_index)
        first_nucleotide_index += 1
    return length


first_nucleotide_index, last_nucleotide_index = 0, (len(sequences[0]) + 1)
sequenceLength(first_nucleotide_index, last_nucleotide_index)

# Create nested dictionary containing frequency of each nucleotide in set of sequences
nucleotide_frequency = {'A': dict.fromkeys(length, 0), 'C': dict.fromkeys(length, 0), 'G': dict.fromkeys(length, 0), 'T': dict.fromkeys(length, 0)}

for sequence in sequences:
    for i, value in enumerate(sequence):
        if sequence[i] == 'A':
            nucleotide_frequency['A'][i] += 1
        elif sequence[i] == 'C':
            nucleotide_frequency['C'][i] += 1
        elif sequence[i] == 'G':
            nucleotide_frequency['G'][i] += 1
        else:
            nucleotide_frequency['T'][i] += 1

# Create array containing nucleotide frequencies
sequence_array = np.zeros(shape=[4, int(len(sequences[0])) + 1])

for i in nucleotide_frequency['A']:
    sequence_array[0][i] = nucleotide_frequency['A'][i]
for i in nucleotide_frequency['C']:
    sequence_array[1][i] = nucleotide_frequency['C'][i]
for i in nucleotide_frequency['G']:
    sequence_array[2][i] = nucleotide_frequency['G'][i]
for i in nucleotide_frequency['T']:
    sequence_array[3][i] = nucleotide_frequency['T'][i]


sequence_array = sequence_array.astype(int)

# Find consensus sequence based on most frequent nucleotide at each position
consensus_sequence_array = sequence_array.argmax(axis=0)

consensus = []
for i, value in enumerate(consensus_sequence_array):
    if consensus_sequence_array[i] == 0:
        consensus.append('A')
    elif consensus_sequence_array[i] == 1:
        consensus.append('C')
    elif consensus_sequence_array[i] == 2:
        consensus.append('G')
    else:
        consensus.append('T')

consensus.pop()
consensus_string = ''.join(consensus)

# Format results for submission
df = pandas.DataFrame(sequence_array)
df.index = ['A:', 'C:', 'G:', 'T:']
df = df.drop(len(sequences[0]), axis = 1)

print(consensus_string)
print(df.to_string(index=True, header=False))
