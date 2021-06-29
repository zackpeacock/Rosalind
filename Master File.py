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
