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