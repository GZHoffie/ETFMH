from db_sketching.utils import Seq2KMers, KMer
from Bio import SeqIO

def compute_resemblence(set_1, set_2):
    intersection = set_1.intersection(set_2)
    union = set_1.union(set_2)
    return len(intersection) / len(union)

def compute_containment(set_1, set_2):
    intersection = set_1.intersection(set_2)
    return len(intersection) / len(set_1)

class KMerSet:
    """
    Base class of the k-mer set.
    Finding the k-mer set given a (list of) sequences.
    Args:
        - k (int): length of the k-mer.
        - canonical (bool): whether we use canonical k-mers.
    """
    def __init__(self, k : int, canonical: bool = False) -> None:
        # The set of k-mers
        self.set = set()

        # Length of k-mer
        self.k = k
        self.canonical = canonical

        # Method to turn sequence into list of k-mers
        self.seq2vec = Seq2KMers(k)

        # Total length of sequences
        self.length = 0
    
    def insert_sequence(self, sequence : str):
        """
        Insert the sequence into k-mer set.
        """
        self.length += len(sequence)
        if self.canonical:
            kmers = self.seq2vec.canonical_kmers(sequence)
        else:
            kmers = self.seq2vec.kmers(sequence)

        self.set.update(kmers)

    def insert_file(self, file):
        """
        Given a fasta file, insert all sequences in the files into the 
        k-mer set.
        Args:
            - file (str): file name of the fasta file .
        """
        for record in SeqIO.parse(file, "fasta"):
            self.insert_sequence(str(record.seq))
    
    def insert_file_list(self, file_list):
        """
        Given a list of fasta files, insert all sequences in the files into the 
        k-mer set.
        Args:
            - file_list (List[str]): a list of fasta files.
        """
        for file in file_list:
            self.insert_file(file)

    def reset(self):
        """
        Set the k-mer set to empty.
        """
        self.set = set()
    
    def resemblence(self, that):
        """
        Calculate the Jaccard index between self.set and that.set.
        Args:
            - that (KMerSet): Another KMerSet object to compare to.
        """
        return compute_resemblence(self.set,that.set)
    
    def containment(self, that):
        """
        Calculate the containment index of self.set in that.set.
        Args:
            - that (KMerSet): Another KMerSet object to compare to.
        """
        return compute_containment(self.set,that.set)

    def ANI_estimation(self, that):
        """
        Estimate ANI using the default estimator.
        Args:
            - that (KMerSet): Another KMerSet object to compare to.
        """
        # Estimator using binomial distribution
        return self.containment(that) ** (1/self.k)
    
        # Estimator using Poisson distribution
        #j = self.containment(that)
        #return 1 + 1/self.k * np.log(2 * j / (1 + j))


class FracMinHash(KMerSet):
    """
    A very simple implementation of FracMinHash.
    Args:
        - condition (function<bool>): a function that takes in a k-mer hash value and
          returns a boolean value. If it returns true, we include that k-mer in our subset.
        - k (int): length of k-mer.
    """
    def __init__(self, condition, k, canonical) -> None:
        super().__init__(k, canonical)
        self.condition = condition
    
    def insert_sequence(self, sequence):
        if self.canonical:
            kmers = self.seq2vec.canonical_kmers(sequence)
        else:
            kmers = self.seq2vec.kmers(sequence)

        self.set.update([i for i in kmers if self.condition(i)])

    

class MeanFracMinHash(KMerSet):
    """
    A version of FracMinHash that stores kmers corresponding to multiple hash function 
    and takes the mean across all of them

    Args:
        - conditons (list[function<bool>]): a list of functions, each of which takes a k-mer hash value and
            returns True if that k-mer should be in the set
        - k (int): length of the k-mer
    """
    def __init__(self, conditions, k, canonical) -> None:
        super().__init__(k, canonical)
        self.conditions = conditions
        self.num_conditions = len(conditions)
        self.sets = [set() for _ in range(self.num_conditions)]
    
    def insert_sequence(self, sequence):
        if self.canonical:
            kmers = self.seq2vec.canonical_kmers(sequence)
        else:
            kmers = self.seq2vec.kmers(sequence)

        for idx, condition in enumerate(self.conditions):
            self.sets[idx].update([i for i in kmers if condition(i)])

    def reset(self):
        """
        Set the k-mer sets to empty.
        """
        self.sets = [set() for _ in range(self.num_conditions)]

    def get_condition_set(self, condition_idx):
        return self.sets[condition_idx]
    
    def resemblence_at_index(self, that, condition_idx):
        """
        Calculate the Jaccard index between self.set and that.set.
        Args:
            - that (MeanFracMinHash): Another MeanFracMinHash object to compare to.
        """
        return compute_resemblence(self.sets[condition_idx],that.sets[condition_idx])
    
    def resemblence(self, that):
        return sum(self.resemblence_at_index(that,idx) for idx in range(self.num_conditions)) / self.num_conditions
    
    def containment_at_index(self, that, condition_idx):
        """
        Calculate the containment index of self.set in that.set.
        Args:
            - that (MeanFracMinHash): Another MeanFracMinHash object to compare to.
        """
        return compute_containment(self.sets[condition_idx],that.sets[condition_idx])

    def containment(self, that):
        return sum(self.containment_at_index(that,idx) for idx in range(self.num_conditions)) / self.num_conditions

    def ANI_estimation(self, that):
        """
        Estimate ANI using the default estimator.
        Args:
            - that (MeanFracMinHash): Another MeanFracMinHash object to compare to.
        """
        # Estimator using binomial distribution
        return self.containment(that) ** (1/self.k)
    
        # Estimator using Poisson distribution
        #j = self.containment(that)
        #return 1 + 1/self.k * np.log(2 * j / (1 + j))
class TruncatedKMerSet(FracMinHash):
    """
    Create (k-l)-mer set based on the constructed k-mer sets and
    try to infer ANI.
    """
    def __init__(self, condition, k, canonical) -> None:
        super().__init__(condition, k, canonical)
    
    def truncate_set(self, l):
        return set([i >> (2 * l) for i in self.set])

    
    def ANI_estimation(self, that):
        # Find containment index of the k-mer set
        kmer_set_containment = len(self.set.intersection(that.set)) / len(self.set)

        # Find containment index of the (k-1)-mer set
        this_k_1_mer_set = set([i >> 2 for i in self.set])
        that_k_1_mer_set = set([i >> 2 for i in that.set])
        k_1_mer_set_containment = len(this_k_1_mer_set.intersection(that_k_1_mer_set)) / len(this_k_1_mer_set)
        print(kmer_set_containment, k_1_mer_set_containment)
        return kmer_set_containment / k_1_mer_set_containment
    
    def k_specific_containment(self, that, k):
        self_low_kmer_set = self.truncate_set(self.k-k)
        that_low_kmer_set = that.truncate_set(that.k-k)
        intersection = len(self_low_kmer_set.intersection(that_low_kmer_set))
        return intersection / len(self_low_kmer_set)

        

class ErrorTolerantFracMinHash(FracMinHash):
    def __init__(self, condition, k, canonical) -> None:
        super().__init__(condition, k, canonical)
        self.kmer = KMer(k)
    
    def insert_sequence(self, sequence):
        if self.canonical:
            kmers = self.seq2vec.canonical_kmers(sequence)
        else:
            kmers = self.seq2vec.kmers(sequence)

        for i in kmers:
            self.set.update([j for j in self.kmer.distance_one_neighbors(i) if self.condition(j)])


if __name__ == "__main__":
    def all(kmer_hash):
        return True
    
    a = TruncatedKMerSet(all, 12)
    a.insert_sequence("CGCGCACGTCGTCGTAC")
    b = TruncatedKMerSet(all, 12)
    b.insert_sequence("CGCGCACGTCGTCGTAG")

    print(a.ANI_estimation(b))