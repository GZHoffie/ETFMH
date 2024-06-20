from db_sketching.kmer_set import KMerSet
from collections import Counter
import numpy as np

class KMerImportanceLearning:
    def __init__(self, kmer_set : KMerSet) -> None:
        self.kmer_set = kmer_set
        self.counter_dict = {}
        self.label_dict = Counter()

    def insert_reference(self, reference_file, label):
        self.kmer_set.reset()
        self.kmer_set.insert_file(reference_file)

        if label not in self.counter_dict:
            self.counter_dict[label] = Counter(set(self.kmer_set.set.keys()))
        else:
            self.counter_dict[label] += Counter(set(self.kmer_set.set.keys()))
        
        self.label_dict[label] += 1

    def _log(self, array):
        """
        Take the log 2 of the elements in the array. Elements that were zero remains zero.
        """
        return np.log2(array, out=np.zeros_like(array, dtype=np.float64), where=(array != 0))
    
    def conditional_entropy(self, kmer: int):
        """
        Find the conditional entropy given a kmer.
        """
        positive_count = []
        negative_count = []

        for label in self.label_dict:
            # Count how many references contain this kmer
            positive_count.append(self.counter_dict[label][kmer])
            negative_count.append(self.label_dict[label] - self.counter_dict[label][kmer])

        positive_count = np.array(positive_count)
        negative_count = np.array(negative_count)

        positive_probs = positive_count / np.sum(positive_count)
        positive_entropy = - np.sum(positive_probs * self._log(positive_probs))

        negative_probs = negative_count / np.sum(negative_count)
        negative_entropy = - np.sum(negative_probs * self._log(negative_probs))


        positive_rate = np.sum(positive_count) / np.sum(positive_count + negative_count)
        res = positive_entropy * positive_rate + negative_entropy * (1-positive_rate)
        return res
    
    def extract_informative_kmers(self, information_gain_threshold = 0):
        res = dict()

        all_kmers = set()
        for label in self.label_dict:
            all_kmers.update(self.counter_dict[label].keys())
        
        orig_probs = np.array(list(self.label_dict.values())) / sum(self.label_dict.values())
        orig_entropy = - np.sum(orig_probs * self._log(orig_probs))

        for kmer in all_kmers:
            information_gain = orig_entropy - self.conditional_entropy(kmer)
            if information_gain > information_gain_threshold:
                res[kmer] = information_gain
        
        return res
    


if __name__ == "__main__":
    from db_sketching.kmer_set import FracMinHash

    