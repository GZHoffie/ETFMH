from db_sketching.utils import Seq2KMers, KMer
from collections import Counter
from Bio import SeqIO
import statistics

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
        - kmer_template (str): a string consisting of 0 and 1's, indicating which base to include in the seed.
              For example, seed_template = "1" * 12 represents a contiguous 12-mer, "1001001" means taking one for every
              3 bases. Or it can be an integer, e.g. `12`, indicating a contiguous 12-mer.
        - canonical (bool): whether we use canonical k-mers.
        - multiplicity (bool): Do we take the multiplicity of k-mers into account when
          calculating the resemblence.
    """
    def __init__(self, kmer_template, canonical: bool = False, multiplicity: bool = False) -> None:
        # The set of k-mers
        self.set = Counter()

        # Length of k-mer
        if isinstance(kmer_template, str):
            self.k = kmer_template.count('1')
        else:
            self.k = kmer_template

        self.canonical = canonical
        self.multiplicity = multiplicity

        # Method to turn sequence into list of k-mers
        if isinstance(kmer_template, str):
            self.seq2vec = Seq2KMers(kmer_template)
        else:
            self.seq2vec = Seq2KMers("1" * kmer_template)

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

    def insert_file(self, file, mode="fasta"):
        """
        Given a fasta file, insert all sequences in the files into the 
        k-mer set.
        Args:
            - file (str): file name of the fasta file .
        """
        for record in SeqIO.parse(file, mode):
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
        self.set = Counter()
    
    def resemblence(self, that):
        """
        Calculate the Jaccard index between self.set and that.set.
        Args:
            - that (KMerSet): Another KMerSet object to compare to.
        """
        if self.multiplicity:
            intersection = sum((self.set & that.set).values())
            union = sum((self.set | that.set).values())
        else:
            intersection = len(self.set & that.set)
            union = len(self.set | that.set)

        return intersection / union
    
    def containment(self, that):
        """
        Calculate the containment index of self.set in that.set.
        Args:
            - that (KMerSet): Another KMerSet object to compare to.
        """
        if self.multiplicity:
            intersection = sum((self.set & that.set).values())
            self_len = sum(self.set.values())
        else:
            intersection = len(self.set & that.set)
            self_len = len(self.set)

        return intersection / self_len

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
    def __init__(self, condition, kmer_template, canonical=False, multiplicity=False) -> None:
        super().__init__(kmer_template, canonical, multiplicity)
        self.condition = condition
    
    def insert_sequence(self, sequence):
        if self.canonical:
            kmers = self.seq2vec.canonical_kmers(sequence)
        else:
            kmers = self.seq2vec.kmers(sequence)

        self.set.update([i for i in kmers if self.condition(i)])

    

class EstimatorFracMinHash(KMerSet):
    """
    A version of FracMinHash that stores kmers corresponding to multiple hash function 
    and uses an estimator to obtain a final result

    Additional Args:
        - conditions (list[function<bool>]): 
            a list of functions, each of which takes a k-mer hash value and
            returns True if that k-mer should be in the set

        - estimator (function from list to float): 
            a function that takes in sampled values and 
            outputs an estimate of the true value
    """
    def __init__(self, conditions, kmer_template, canonical=False, multiplicity=False, estimator=statistics.mean) -> None:
        super().__init__(kmer_template, canonical, multiplicity)
        self.conditions = conditions
        self.num_conditions = len(conditions)
        self.sets = [set() for _ in range(self.num_conditions)]
        self.estimator = estimator
    
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
        return self.estimator(self.resemblence_at_index(that,idx) for idx in range(self.num_conditions))
    
    def containment_at_index(self, that, condition_idx):
        """
        Calculate the containment index of self.set in that.set.
        Args:
            - that (MeanFracMinHash): Another MeanFracMinHash object to compare to.
        """
        return compute_containment(self.sets[condition_idx],that.sets[condition_idx])

    def containment(self, that):
        return self.estimator(self.containment_at_index(that,idx) for idx in range(self.num_conditions))

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


class MultiSeedFracMinHash(FracMinHash):
    pass

class TruncatedKMerSet(FracMinHash):
    """
    Create (k-l)-mer set based on the constructed k-mer sets and
    try to infer ANI.
    """
    def __init__(self, condition, kmer_template) -> None:
        super().__init__(condition, kmer_template, canonical=False, multiplicity=False)
        
    
    def insert_sequence(self, sequence):
        # Store only the original kmers
        kmers = self.seq2vec.kmers(sequence)

        self.set.update([i for i in kmers if self.condition(i)])
    
    def truncate_set(self, l):
        return set([i >> (2 * l) for i in self.set])
    
    def conditional_probability(self, that):
        """
        Calculate Pr[this base is not mutated | previous k-1 bases are not mutated] using
        the truncated k-mer set.
        """
        # Find containment index of the (k-1)-mer set
        this_k_1_mer_set = set([i >> 2 for i in self.set.keys()])
        that_k_1_mer_set = set([i >> 2 for i in that.set.keys()])

        # Find containment index of the k-mer set
        kmer_set_containment = len(self.set & that.set) / len(self.set)
        k_1_mer_set_containment = len(this_k_1_mer_set.intersection(that_k_1_mer_set)) / len(this_k_1_mer_set)
        return kmer_set_containment / k_1_mer_set_containment
    
    
    def ANI_estimation(self, that):
        # Find containment index of the (k-1)-mer set
        this_k_1_mer_set = set([i >> 2 for i in self.set.keys()])
        that_k_1_mer_set = set([i >> 2 for i in that.set.keys()])

        # Find containment index of the k-mer set
        kmer_set_containment = len(self.set & that.set) / len(self.set)
        k_1_mer_set_containment = len(this_k_1_mer_set.intersection(that_k_1_mer_set)) / len(this_k_1_mer_set)
        print("k-mer containment", kmer_set_containment, "(k-1)-mer set containment", k_1_mer_set_containment)
        p1 = kmer_set_containment ** (1/(self.k))
        p2 = k_1_mer_set_containment ** (1/(self.k-1))
        p3 = kmer_set_containment / k_1_mer_set_containment
        print("Estimation using k", kmer_set_containment ** (1/(self.k)))
        print("Estimation using k-1", k_1_mer_set_containment ** (1/(self.k-1)))
        print("Estimation using conditional prob", kmer_set_containment / k_1_mer_set_containment)
        return (p1 + p2 + p3) / 3
    
    def k_specific_containment(self, that, k):
        self_low_kmer_set = self.truncate_set(self.k-k)
        that_low_kmer_set = that.truncate_set(that.k-k)
        intersection = len(self_low_kmer_set.intersection(that_low_kmer_set))
        return intersection / len(self_low_kmer_set)

        

class ErrorTolerantFracMinHash(FracMinHash):
    def __init__(self, condition, kmer_template, canonical=False) -> None:
        super().__init__(condition, kmer_template, canonical, multiplicity=False)
        self.kmer = KMer(kmer_template)
    
    def insert_sequence(self, sequence):
        if self.canonical:
            kmers = self.seq2vec.canonical_kmers(sequence)
        else:
            kmers = self.seq2vec.kmers(sequence)

        for i in kmers:
            self.set.update([j for j in self.kmer.distance_one_neighbors(i) if self.condition(j)])


class NullomerSet(KMerSet):
    def __init__(self, kmer_template, canonical: bool = False, multiplicity: bool = False) -> None:
        super().__init__(kmer_template, canonical, multiplicity)
        self.set = self._init_set()
    

    def _init_set(self):
        canonical_kmer_dict = self.seq2vec.all_canonical_kmers()
        res = set(canonical_kmer_dict.values())
        return res
    
    def insert_sequence(self, sequence : str):
        """
        Insert the sequence into k-mer set.
        """
        self.length += len(sequence)
        if self.canonical:
            kmers = self.seq2vec.canonical_kmers(sequence)
        else:
            kmers = self.seq2vec.kmers(sequence)

        self.set.remove(kmers)
        


if __name__ == "__main__":
    def all(kmer_hash):
        return True
    
    a = KMerSet(31, True, True)
    a.insert_file("/home/zhenhao/tc-benchmark/data/sensitivity_test/s_aureus_coverage_0.0535.fastq", "fastq")
    N3 = len([i for i in a.set if a.set[i] == 3])
    N2 = len([i for i in a.set if a.set[i] == 2])
    N1 = len([i for i in a.set if a.set[i] == 1])

    print(N3, N2, N1, N2/N1 * 2, N3/N2 * 3)