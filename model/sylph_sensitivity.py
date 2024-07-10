import numpy as np

true_coverage_list = np.arange(0.01, 0.02, 0.01)
read_error_rate = 0.05
num_kmers = 4000000 # for E. coli
num_simulations = 100
k = 21

for coverage in true_coverage_list:
    # Find effective coverage
    effective_coverage = coverage * ((1-read_error_rate) ** k)
    num_reported = 0
    total_num = 0
    for _ in range(num_simulations):
        samples = np.random.poisson(effective_coverage, num_kmers)
        containment_index = np.sum(samples > 0) / len(samples)
        N1 = np.sum(samples == 1)
        N2 = np.sum(samples == 2)
        if N1 < 3 or N2 < 3:
            continue
        estimated_lambda = N2 / N1
        estimated_ANI = (containment_index / (1-np.exp(-estimated_lambda))) ** (1/k)
        print(N2, N1, containment_index, estimated_ANI)
        if estimated_ANI > 0.95:
            num_reported += 1
        total_num += 1

    print(effective_coverage, num_reported / total_num, total_num)



