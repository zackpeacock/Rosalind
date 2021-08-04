import itertools

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


seqs_dict = FormatFASTA('rosalind_long.txt')
seqs_list = list(seqs_dict.values())
seqs = list(itertools.permutations(seqs_list))


def glue(seqs):
    sequences = []
    for i in range(len(seqs)):
        sequence = seqs[i][0]
        for j in range(1, len(seqs[i])):
            next_seq = seqs[i][j]
            m = int((len(next_seq) / 2))
            for k in range(m + 1, len(next_seq)):
                if next_seq[0:k] == sequence[-k:]:
                    overlap = next_seq[0:k]
                    seq = next_seq.replace(overlap, '')
                    sequence = sequence + seq
                    sequences.append(sequence)

    sequences.sort(key=len, reverse=True)
    superstring = sequences[0]
    return superstring