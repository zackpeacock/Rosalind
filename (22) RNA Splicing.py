# Given: A DNA string s (of length at most 1 kbp) and a collection of substrings of s acting as introns. All strings are given in FASTA format.

# Return: A protein string resulting from transcribing and translating the exons of s. 


# Function to create dictionary of DNA sequence and introns with accompanying FASTA IDs
def FormatFASTA(filename):
    file = open(filename, 'r')
    sequences_file = [line.strip() for line in file.readlines()]

    sequences = {}
    for line in sequences_file:
        if '>' in line:
            sequence_label = line
            sequences[sequence_label] = ''
        else:
            sequences[sequence_label] += line
    return sequences


sequences_dict = FormatFASTA('rosalind_splc.txt')

# Identify DNA sequence, assign to variable, then remove from dictionary
for i in sequences_dict:
    if len(sequences_dict[i]) == max(len(sequence) for sequence in sequences_dict.values()):
        DNA_sequence = sequences_dict[i]
        del sequences_dict[i]
        break

# Sequences dictionary now exclusively contains introns
introns = sequences_dict

# Copy/paste previously formulated dictionary of RNA codons with amino acid counterparts
RNA_codons = {
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
    'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'
}

# Create "DNA codon" dictionary by replacing uracil with thymine
DNA_codons = {codon.replace('U', 'T'): aa for codon, aa in RNA_codons.items()}

# Splice introns from DNA sequence
spliced_sequence = DNA_sequence
for i in introns:
    if introns[i] in DNA_sequence:
        spliced_sequence = spliced_sequence.replace(introns[i], '')

# Translate DNA sequence
protein_string = ''
for i in range(0, len(spliced_sequence) - 3, 3):
    codon = spliced_sequence[i:i + 3]
    protein_string = protein_string + str(DNA_codons[codon])

print(protein_string)
