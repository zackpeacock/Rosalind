# Given: A positive integer n ≤ 10000 followed by a permutation π of length n.

# Return: A longest increasing subsequence of π, followed by a longest decreasing subsequence of π.


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
