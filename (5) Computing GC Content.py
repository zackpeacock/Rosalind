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