# Given: Positive integers n ≤ 100 and m ≤ 20.

# Return: The total number of pairs of rabbits that will remain after the n-th month if all rabbits live for m months.


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
