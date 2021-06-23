# Given: A DNA string of length at most 1000 nt.

# Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in the DNA string.


nucleotides = {}
DNAsequence = open('rosalind_dna.txt', 'r')

# Create dictionary of nucleotides and corresponding frequencies
for sequence in DNAsequence:
    for nucleotide in sequence:
        if nucleotide not in nucleotides:
            nucleotides[nucleotide] = 1
        else:
            nucleotides[nucleotide] += 1
print(str(nucleotides['A']) + ' ' + str(nucleotides['C']) + ' ' + str(nucleotides['G']) + ' ' + str(nucleotides['T']))
