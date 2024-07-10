from db_sketching.kmer_set import KMerSet
import random

class GenomeClassifer:
    def __init__(self, kmer_set) -> None:
        self.kmer_set = kmer_set

    def predict(self, genome_file_1, genome_file_2):
        kmer_set_1 = self.kmer_set()
        kmer_set_2 = self.kmer_set()

        # Read the genome file
        kmer_set_1.insert_file(genome_file_1)
        kmer_set_2.insert_file(genome_file_2)

        return kmer_set_1.resemblence(kmer_set_2)
    
    def AUROC(self, genome_file_dict, sample_pairs=100, positive_rate=0.5):
        # Set of all taxonomy labels
        label_set = set(genome_file_dict)


        for _ in range(sample_pairs):
            if random.random() > positive_rate:
                # Sample pair from the same 

