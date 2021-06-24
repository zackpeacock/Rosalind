# Given: An RNA string corresponding to a strand of mRNA (of length at most 10 kbp).

# Return: The encoded protein string.


file = open('rosalind_prot.txt', 'r')
RNAsequence = file.read()

# Create list of codons contained within RNA sequence
codons = []
n = 3
for i in range(0, len(RNAsequence), n):
    codons.append(RNAsequence[i:i+n])

# Create dictionary of codons with corresponding amino acids
codon_aminoAcid_pairs = {
'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V',
'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V',
'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
'UAA': 'Stop', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
'UAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
'UGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G',
}

# Convert codons to amino acids
aminoAcids = []
for codon in codons:
    if codon in codon_aminoAcid_pairs:
        aminoAcids.append(codon_aminoAcid_pairs[codon])

# Convert amino acid list to polypeptide string and remove 'Stop'
polypeptide = ''.join(aminoAcids)
if polypeptide[-4:len(polypeptide)] == 'Stop':
    polypeptide = polypeptide.replace(polypeptide[-4:len(polypeptide)], '')

print(polypeptide)
