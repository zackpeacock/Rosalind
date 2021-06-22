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
