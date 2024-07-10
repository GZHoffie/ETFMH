import numpy as np
from collections import Counter

true_coverage_list = np.arange(0, 0.3, 0.005)#[0.008, 0.015, 0.0284, 0.0535, 0.101, 0.1902, 0.3585, 0.6757]
read_error_rate = 0.05
num_kmers = 4000000 # for E. coli
c = 200
num_simulations = 1000
k = 21

for coverage in true_coverage_list:
    # Find effective coverage
    effective_coverage = coverage #* ((1-read_error_rate) ** k)
    num_reported = 0
    total_num = 0
    for _ in range(num_simulations):
        total_num += 1

        # Sample a number of k-mers
        kmer_samples = np.random.choice(np.arange(int(num_kmers / c)), int(effective_coverage * num_kmers / c))
        kmer_counter = Counter(kmer_samples)
        
        # Find the containment index
        containment_index = len(kmer_counter) / int(num_kmers / c)
        N1 = list(kmer_counter.values()).count(1)
        N2 = list(kmer_counter.values()).count(2)
        if N1 < 3 or N2 < 3:
            estimated_lambda = 1
        else:
            estimated_lambda = N2 / N1
        estimated_ANI = (containment_index / (1-np.exp(-estimated_lambda))) ** (1/k)
        #print(effective_coverage, containment_index, estimated_lambda, estimated_ANI)
        if estimated_ANI > 0.95:
            num_reported += 1
        

    print(effective_coverage, num_reported / total_num, total_num)



