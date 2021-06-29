# Given: At most 15 UniProt Protein Database access IDs.

# Return: For each protein possessing the N-glycosylation motif ( N{P}[ST]{P} ), output its given access ID followed by a list of locations in the protein string
# where the motif can be found.


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
