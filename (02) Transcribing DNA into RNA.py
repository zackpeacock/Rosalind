# Given: A DNA string having length at most 1000 nt.

# Return: The transcribed RNA string of the DNA string.


DNAsequence = open('rosalind_rna.txt', 'r')

# Substitute uracil for thymine to transcribe sequence
for sequence in DNAsequence:
    RNAsequence = sequence.replace('T', 'U')
print(RNAsequence)
