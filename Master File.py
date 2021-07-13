# *** Counting DNA Nucleotides (DNA) ***

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


# *** Transcribing DNA into RNA (RNA) ***

DNAsequence = open('rosalind_rna.txt', 'r')

# Substitute uracil for thymine to transcribe sequence
for sequence in DNAsequence:
    RNAsequence = sequence.replace('T', 'U')
print(RNAsequence)


# *** Complementing a Strand of DNA (REVC) ***

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


# *** Rabbits and Recurrence Relations (FIB) ***

file = open('rosalind_fib.txt', 'r')

# Assign months of mating and offspring per month to variables
givenIntegers = file.read().strip().split(' ')
month = int(givenIntegers[0])
offspring = int(givenIntegers[1])

# Calculate final number of rabbit pairs
def rabbitPairs(month, offspring):
    if month <= 2:
        return 1
    else:
        return rabbitPairs(month - 1, offspring) + (offspring * (rabbitPairs(month - 2, offspring)))


print(rabbitPairs(month, offspring))


# *** Computing GC Content (GC) ***

# Remove newlines from FASTA file
def formatFASTA(filename):
    file = open(filename, 'r')
    return [line.strip() for line in file.readlines()]

# Create dictionary of FASTA labels with accompanying sequences
sequencesFile = formatFASTA('rosalind_gc.txt')
sequenceLabel = ''
sequencesDict = {}

for line in sequencesFile:
    if '>' in line:
        sequenceLabel = line
        sequencesDict[sequenceLabel] = ''
    else:
        sequencesDict[sequenceLabel] += line

# Computes % GC content of sequence
def compute_GC_content(sequence):
    return round(((sequence.count('G') + sequence.count('C')) / len(sequence) * 100), 5)

# Create dictionary of FASTA labels with accompanying GC content
GC_content_dict = {}
for label, sequence in sequencesDict.items():
    GC_content_dict[label] = compute_GC_content(sequence)

# Find FASTA label with highest GC content and return label/content pair
max_GC_label = max(GC_content_dict, key = GC_content_dict.get)
max_GC_content = max_GC_label.replace('>', '') + '\n' + str(GC_content_dict[max_GC_label])

print(max_GC_content)


# *** Counting Point Mutations (HAMM) ***

file = open('rosalind_hamm.txt', 'r')

# Assign sequences to be compared to variables and initialize mutation counter
rawSequences = file.read().split('\n')
originalSequence = rawSequences[0]
mutatedSequence = rawSequences[1]
mutations = 0

# Compute number of mutations observed between sequences
for i, value in enumerate(originalSequence):
    if originalSequence[i] == mutatedSequence[i]:
        pass
    else:
        mutations += 1

print(mutations)


# *** Mendel's First Law (IPRB) ***

integers = open('rosalind_iprb.txt', 'r')
for string in integers:
    integerList = string.strip().split(' ')

# Assign numbers of homozygous dominant (k), heterozygous (m), and homozygous recessive (n) individuals to variables
k = int(integerList[0])
m = int(integerList[1])
n = int(integerList[2])

# t = population total
t = k + m + n

# Probabilities k, m, or n is first selection
probk1 = k / t
probm1 = m / t
probn1 = n / t

# Probabilities k, m, or n is second selection if not first selection
probk2 = k / (t - 1)
probm2 = m / (t - 1)
probn2 = n / (t - 1)

# Probabilities k, m, or n is second selection if also first selection
probk3 = (k - 1) / (t - 1)
probm3 = (m - 1) / (t - 1)
probn3 = (n - 1) / (t - 1)

# Probability of two recessive alleles
probRec = (0.25 * (probm1 * probm3)) + (0.5 * (probm1 * probn2)) + (probn1 * probn3) + (0.5 * (probn1 * probm2))

# Probability of one dominant allele
probDom = 1 - probRec

print(round(probDom, 5))


# *** Translating RNA Into Protein (PROT) ***

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


# *** Finding a Motif in DNA (SUBS) ***

file = open('rosalind_subs.txt', 'r')

# Assign sequence and subsequence to variables
rawSequences = file.read().split('\n')
sequence = rawSequences[0]
subsequence = rawSequences[1]

# Create list of locations of subsequence
locations = []
n = len(subsequence)
for i, value in enumerate(sequence):
    if sequence[i:i+n] == subsequence:
        locations.append(str(i+1))

print(' '.join(locations))


# *** Consensus and Profile (CONS) ***
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


# *** Mortal Fibonacci Rabbits (FIBD) ***

file = open('rosalind_fibd.txt', 'r')

# Assign total number of months and lifespan of rabbits to variables
givenIntegers = file.read().strip().split(' ')
total_months = int(givenIntegers[0])
lifespan = int(givenIntegers[1])

# Create list of rabbit pairs in each month
pairs = [1, 1]
month = 2
while month < total_months:

    # Normal fibonacci progression prior to any rabbits dying
    if month < lifespan:
        pairs.append(pairs[-1] + pairs[-2])

    # Recursively calculate rabbit pairs per month once rabbits begin to die
    elif month == lifespan or month == lifespan + 1:
        pairs.append(pairs[-1] + pairs[-2] - 1)
    else:
        pairs.append(pairs[-1] + pairs[-2] - pairs[month - (lifespan + 1)])
    month += 1

print(pairs[-1])


# *** Overlap Graphs (GRPH) ***

# Remove newlines from FASTA file
def formatFASTA(filename):
    file = open(filename, 'r')
    return [line.strip() for line in file.readlines()]

# Create dictionary of FASTA labels with accompanying sequences
sequencesFile = formatFASTA('rosalind_grph.txt')
sequenceLabel = ''
sequencesDict = {}

for line in sequencesFile:
    if '>' in line:
        sequenceLabel = line
        sequencesDict[sequenceLabel] = ''
    else:
        sequencesDict[sequenceLabel] += line

# Create suffix and prefix dictionaries for each sequence, matched with its SeqID
suffixes = {}
prefixes = {}

for key in sequencesDict:
    suffixes[key] = sequencesDict[key][-3:]
    prefixes[key] = sequencesDict[key][:3]

# Identify suffixes that match prefixes and sort SeqIDs into separate lists
suffix_ids = []
prefix_ids = []
for key in suffixes:
    for key2 in prefixes:
        if suffixes[key] == prefixes[key2]:
            if key != key2:
                suffix_ids.append(key.replace('>', ''))
                prefix_ids.append(key2.replace('>', ''))

# Return SeqID pairs for matches
for i in range(len(suffix_ids)):
    print(suffix_ids[i] + ' ' + prefix_ids[i])


# *** Calculating Expected Offspring (IEV) ***

# Store given values in list
file = open('rosalind_iev.txt', 'r')
givenIntegers = file.read().strip().split(' ')

# Create dictionary with probabilities that each type of couple produces offspring with the dominant phenotype
prob_dom_phenotype = {'AA-AA': 1, 'AA-Aa': 1, 'AA-aa': 1, 'Aa-Aa': 0.75, 'Aa-aa': 0.5, 'aa-aa': 0}

# Create dictionary with number of each type of couple
num_couples = {'AA-AA': int(givenIntegers[0]), 'AA-Aa': int(givenIntegers[1]), 'AA-aa': int(givenIntegers[2]),
               'Aa-Aa': int(givenIntegers[3]), 'Aa-aa': int(givenIntegers[4]), 'aa-aa': int(givenIntegers[5])}

# Multiply probability of dominant phenotype by number of each type of couple and store in dictionary
num_dom_phenotype = {}
for key in prob_dom_phenotype:
    for key2 in num_couples:
        if key == key2:
            # Multiply by 2 to account for 2 offspring per couple
            num_dom_phenotype[key] = 2 * prob_dom_phenotype[key] * num_couples[key2]

# Sum number of offspring with dominant phenotype expected to be produced by each type of couple
total_dom_phenotype = sum(num_dom_phenotype.values())

print(total_dom_phenotype)


# *** Finding a Shared Motif (LCSM) ***

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


# *** Independent Alleles (LIA) ***

# Assign given values for k and N to variables
import math
file = open('rosalind_lia.txt', 'r')
givenIntegers = file.read().strip().split(' ')

k = int(givenIntegers[0])
N = int(givenIntegers[1])
P = 2 ** k

# Probability of heterozygous offspring if at least 1 parent is heterozygous = 0.25
# Probability of N heterozygotes in total population P = [P! / N!(P - N)!] * (0.25 ^ i) * [0.75 ^ (P - i)]
# Sum probabilities of producing N, N + 1, N + 2 ... P heterozygotes
probability = 0
for i in range(N, P + 1):
    single_gen_prob = (math.factorial(P) / (math.factorial(i) * math.factorial(P - i))) * (0.25 ** i) * (0.75 ** (P - i))
    probability += single_gen_prob

print(round(probability, 3))


# *** Finding a Protein Motif (MPRT) ***

import requests

# Compile given protein IDs into list
ID_list = open('rosalind_mprt.txt', 'r')
proteinIDs = ID_list.read().split('\n')

# Write FASTA text into .txt file
proteins = dict.fromkeys(proteinIDs)
for ID in proteins.keys():
    protein_fasta = requests.get('https://uniprot.org/uniprot/' + str(ID) + '.fasta')
    with open(str(ID) + '.txt', 'w') as file:
        file.write(protein_fasta.text)
        file.close()

# Create dictionary of protein IDs paired with their respective amino acid sequences
def formatFASTA(filename):
    file = open(filename, 'r')
    return [line.strip() for line in file.readlines()]


for ID in proteins.keys():
    fragments = formatFASTA(ID + '.txt')
    fragments.remove(fragments[0])
    sequence = ''.join(fragments)
    proteins[ID] = sequence

# Create dictionary of protein IDs paired with indices at which N-glycosylation motif is found
motif_occurrences = {}
for ID in proteins.keys():
    motif_occurrences.setdefault(ID, [])
    seq = proteins[ID]
    for i in range(len(seq) - 3):
        if seq[i] == 'N':
            if seq[i + 1] != 'P':
                if seq[i + 2] == 'S':
                    if seq[i + 3] != 'P':
                        motif_occurrences[ID].append(str(i + 1))
                elif seq[i + 2] == 'T':
                    if seq[i + 3] != 'P':
                        motif_occurrences[ID].append(str(i + 1))

# Only print protein IDs that have N-glycosylation motif
for ID in motif_occurrences.keys():
    if motif_occurrences[ID] != []:
        print(ID + '\n' + ' '.join(motif_occurrences[ID]))


# *** Inferring mRNA from Protein (MRNA) ***

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


# *** Open Reading Frames (ORF) ***

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


# *** Enumerating Gene Orders (PERM) ***

from itertools import permutations
import math

# Assign given integer to variable n
file = open('rosalind_perm.txt', 'r')
given_integer = file.read().split()
n = int(given_integer[0])

# Create list of permutations of n
perms = list(permutations(range(1, n + 1)))

# Format list of permutations
formatted_perms = []
for lst in perms:
    sublist = []
    for i in lst:
        sublist.append(str(i))
    formatted_perms.append(sublist)

# Print number of permutations and all possible permutations in proper format
print(math.factorial(n))
for lst in formatted_perms:
    print(' '.join(lst))


# *** Calculating Protein Mass (PRTM) ***

file = open('rosalind_prtm.txt', 'r')
protein_sequence = file.read().strip()

# Create dictionary of amino acid masses from table provided at rosalind.info
amino_acid_masses = {
    'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841, 'G': 57.02146, 'H': 137.05891,
    'I': 113.08406, 'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293, 'P': 97.05276, 'Q': 128.05858,
    'R': 156.10111, 'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333
}

# Sum weights of amino acids in protein sequence
weight = 0
for amino_acid in protein_sequence:
    weight += amino_acid_masses[amino_acid]

# Round weight to 3 decimal places
print(round(weight, 3))


# *** Locating Restriction Sites (REVP) ***

# Function to isolate DNA sequence from FASTA file
def IsolateSequence(filename):
    file = open(filename, 'r')
    sequence_file = [line.strip() for line in file.readlines()]
    DNA_sequence = ''
    for line in sequence_file:
        if '>' in line:
            pass
        else:
            DNA_sequence += line
    return DNA_sequence

# Function to produce reverse complement of DNA sequence
def ReverseComplement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complement = ''
    for nucleotide in sequence:
        reverse_complement = complement[nucleotide] + reverse_complement
    return reverse_complement

# Assign provided sequence to variable
DNA_sequence = IsolateSequence('rosalind_revp.txt')

# Create dictionary of locations where reverse complement of DNA substring is the same as the substring itself
matches = {}
for i in range(len(DNA_sequence) - 3):
    for j in range(i + 4, i + 13):
        if j > len(DNA_sequence):
            break
        else:
            substring = DNA_sequence[i:j]
            if substring == ReverseComplement(substring):
                matches[i + 1] = substring

# Print locations and lengths of substrings
for i in matches:
    print(str(i) + ' ' + str(len(matches[i])))


# *** RNA Splicing (SPLC) ***

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


# *** Enumerating k-mers Lexicographically ***

import itertools

# Create and format alphabet list from .txt file, assign string length to variable
file = open('rosalind_lexf.txt', 'r')
alphabet = file.read().strip().split(' ')
length = alphabet[-1][-1]
alphabet[-1] = alphabet[-1][0]

# Create list of all possible permutations of given string length, including those in which letters repeat
# Python automatically sorts lexicographically
permutations = []
for i in itertools.product(alphabet, repeat=int(length)):
    permutations.append(i)

# Convert permutations to string form
strings = []
for i in permutations:
    string = ''.join(i)
    strings.append(string)

# Print individual strings
for string in strings:
    print(string)


# *** Longest Increasing Subsequence (LGIS) ***

# Format permutation into list
file = open('rosalind_lgis.txt', 'r')
raw_permutation = file.read().strip().split(' ')
raw_permutation[0] = raw_permutation[0][5:]

# Convert values from int to string
permutation = []
for i in raw_permutation:
    permutation.append(int(i))


# Function to determine longest increasing subsequence of given permutation
def LongestIncreasingSubsequence(sequence):
    l = len(sequence)

    # Initialize list to populate with second to last index
    prev_index = [0] * l
    for i in range(l):
        prev_index[i] = i

    # Initialize list to be populated with LIS candidates
    subsequence = [1] * l

    # Find LIS candidates
    for i in range(1, l):
        for j in range(i):
            if sequence[i] > sequence[j] and subsequence[i] < subsequence[j] + 1:
                subsequence[i] = subsequence[j] + 1
                prev_index[i] = j

    # Initialize maximum LIS candidate and index of maximum LIS candidate
    max = 0
    max_index = 0

    # Find maximum LIS candidate
    for i in range(l):
        if max < subsequence[i]:
            max = subsequence[i]
            max_index = i

    LIS = [sequence[max_index]]
    while max_index != prev_index[max_index]:
        max_index = prev_index[max_index]
        LIS.append(sequence[max_index])

    return max, reversed(LIS)


# Function to determine longest decreasing subsequence of given permutation
def LongestDecreasingSubsequence(sequence):
    l = len(sequence)

    # Initialize list to populate with second to last index
    prev_index = [0] * l
    for i in range(l):
        prev_index[i] = i

    # Initialize list to be populated with LDS candidates
    subsequence = [1] * l

    # Find LDS candidates
    for i in range(1, l):
        for j in range(i):
            if sequence[i] < sequence[j] and subsequence[i] < subsequence[j] + 1:
                subsequence[i] = subsequence[j] + 1
                prev_index[i] = j

    # Initialize maximum LDS candidate and index of maximum LDS candidate
    max = 0
    max_index = 0

    # Find maximum LDS candidate
    for i in range(l):
        if max < subsequence[i]:
            max = subsequence[i]
            max_index = i

    LDS = [sequence[max_index]]
    while max_index != prev_index[max_index]:
        max_index = prev_index[max_index]
        LDS.append(sequence[max_index])

    return max, reversed(LDS)


LIS = LongestIncreasingSubsequence(permutation)
LDS = LongestDecreasingSubsequence(permutation)

# Print LIS and LDS formatted for submission
print(' '.join(str(i) for i in LIS[1]))
print(' '.join(str(i) for i in LDS[1]))