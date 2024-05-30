from db_sketching.utils import Seq2KMers
from Bio import SeqIO
import numpy as np

class KMerSet:
    def __init__(self, k) -> None:
        self.set = set()
        self.k = k
        self.seq2vec = Seq2KMers(k)
    
    def insert_sequence(self, sequence):
        kmers = self.seq2vec.canonical_kmers(sequence)
        self.set.update(kmers)
    
    def insert_file_list(self, file_list):
        for file in file_list:
            for record in SeqIO.parse(file, "fasta"):
                self.insert_sequence(str(record.seq))

    def reset(self):
        self.set = set()
    
    def resemblence(self, that):
        intersection = len(self.set.intersection(that.set))
        union = len(self.set.union(that.set))
        return intersection / union
    
    def containment(self, that):
        intersection = len(self.set.intersection(that.set))
        return intersection / len(self.set)

    def ANI_estimation(self, that):
        #return self.containment(that) ** (1/self.k)
        j = self.resemblence(that)
        return 1 + 1/self.k * np.log(2 * j / (1 + j))


class FracMinHash(KMerSet):
    def __init__(self, condition, k) -> None:
        super().__init__(k)
        self.condition = condition
    
    def insert_sequence(self, sequence):
        kmers = self.seq2vec.canonical_kmers(sequence)
        self.set.update([i for i in kmers if self.condition(i)])


class ErrorTolerantFracMinHash(FracMinHash):
    def __init__(self, condition, k) -> None:
        super().__init__(condition, k)
        self.error_tolerance_map = self._init_map()

