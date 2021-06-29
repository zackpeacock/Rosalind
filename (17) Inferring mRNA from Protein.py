file = open('rosalind_mrna.txt', 'r')
protein = file.read().strip()

# Copy/paste preexisting codon: amino acid dictionary
codons = {
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

# Invert codon dictionary such that keys are amino acids and values are codons
inverted_codons = {}
for i in codons:
    for j in codons:
        if codons[i] == codons[j]:
            inverted_codons.setdefault(codons[i], [i]).append(j)

# Remove duplicate codons
for i in inverted_codons:
    inverted_codons[i] = list(set(inverted_codons[i]))

# Multiply number of codons for each amino acid
possible_strings = 1
for aa in protein:
    possible_strings = possible_strings * (len(inverted_codons[aa]))

# Finally, multiply by 3 to account for 3 stop codons
print((possible_strings * 3) % 1000000)
