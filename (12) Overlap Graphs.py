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
