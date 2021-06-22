DNAsequence = open('rosalind_rna.txt', 'r')

# Substitute uracil for thymine to transcribe sequence
for sequence in DNAsequence:
    RNAsequence = sequence.replace('T', 'U')
print(RNAsequence)
