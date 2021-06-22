DNAsequence = open('rosalind_revc.txt', 'r')

# Create dictionary of complementary base pairs
complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# Find reverse complement of sequence
reverseComplement = ''
for sequence in DNAsequence:
    for nucleotide in sequence:
        if nucleotide == '\n':
            pass
        else:
            reverseComplement = complement[nucleotide] + reverseComplement
print(reverseComplement)
