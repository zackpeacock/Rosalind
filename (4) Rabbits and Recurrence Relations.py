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