# Isolate DNA sequence from FASTA file
def formatFASTA(filename):
    file = open(filename, 'r')
    return [line.strip() for line in file.readlines()]


sequence_file = formatFASTA('rosalind_orf.txt')
DNA_sequence = ''
for line in sequence_file:
    if '>' in line:
        pass
    else:
        DNA_sequence += line

# Copy/paste preexisting codon: amino acid dictionary
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

# Define function to determine viable protein strings resulting from all possible reading frames of a sequence
def OpenReadingFrames(sequence):

    # Find indices of start codons
    start_indices = []
    for i in range(len(sequence) - 2):
        if sequence[i:i + 3] == 'ATG':
            start_indices.append(i)

    # Create dictionary of strings corresponding to each reading frame
    strings = {}
    for i in start_indices:
        strings.setdefault(i, [])
        for j in range(i, len(sequence) - 2, 3):
            if DNA_codons[sequence[j:j + 3]] == 'Stop':
                strings[i].append(DNA_codons[sequence[j:j + 3]])
                break
            else:
                strings[i].append(DNA_codons[sequence[j:j + 3]])

    for key in strings:
        strings[key] = ''.join(strings[key])

    # Remove strings that do not end with a stop codon
    viable_strings = []
    for value in strings.values():
        if value[-4:] == 'Stop':
            value = value[:-4]
            viable_strings.append(value)
    return viable_strings


# Find reverse complement of DNA sequence
complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
reverse_complement = ''
for i in DNA_sequence:
    reverse_complement = complement[i] + reverse_complement

# Create list of all possible protein candidates resulting from DNA sequence and its reverse complement
protein_candidates = []
for string in OpenReadingFrames(DNA_sequence):
    protein_candidates.append(string)
for string in OpenReadingFrames(reverse_complement):
    protein_candidates.append(string)

# Remove duplicates from list
edited_protein_candidates = list(set(protein_candidates))

# Print candidates
for candidate in edited_protein_candidates:
    print(candidate)
