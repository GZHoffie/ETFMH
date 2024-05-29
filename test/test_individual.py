from db_sketching.kmer_set import KMerSet, FracMinHash

data_1 = ["./data/562.fna"]
data_2 = ["./data/564.fna"]

k = 31

full1 = KMerSet(k)
full2 = KMerSet(k)
def cond(kmer_hash):
    hash = (976369 * kmer_hash + 1982627) % 10000
    if hash < 100:
        return True
    else:
        return False
frac1 = FracMinHash(cond, k)
frac2 = FracMinHash(cond, k)

# insert dataset
full1.insert_file_list(data_1)
full2.insert_file_list(data_2)
frac1.insert_file_list(data_1)
frac2.insert_file_list(data_2)

print(full1.resemblence(full2))
print(frac1.resemblence(frac2))

print(full1.ANI_estimation(full2))
print(frac1.ANI_estimation(frac2))

