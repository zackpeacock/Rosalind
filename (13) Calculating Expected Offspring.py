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