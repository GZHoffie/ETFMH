from db_sketching.kmer_set import KMerSet
from collections import Counter
import numpy as np
import pickle

class KMerSelection:
    def __init__(self):
        pass

    def select_kmers(self, **kwargs):
        """
        Return a dict of k-mers that are helpful based on the scheme.
        """
        pass

class KMerInformationGain(KMerSelection):
    """
    Select the set of k-mers that have high information gain.
    """
    def __init__(self, kmer_set : KMerSet) -> None:
        super().__init__()
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

    def store(self, file_name):
        """
        Store self.counter_dict and self.label_dict in a pickle file "file_name".
        """
        import pickle

        with open(file_name, "wb") as f:
            pickle.dump((self.counter_dict, self.label_dict), f)
    
    def load(self, file_name):
        """
        Load the pickle file into self.counter_dict and self.label_dict.
        """
        import pickle

        with open(file_name, "rb") as f:
            self.counter_dict, self.label_dict = pickle.load(f)


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
    
    def select_kmers(self, **kwargs):
        res = dict()

        all_kmers = set()
        for label in self.label_dict:
            all_kmers.update(self.counter_dict[label].keys())
        
        orig_probs = np.array(list(self.label_dict.values())) / sum(self.label_dict.values())
        orig_entropy = - np.sum(orig_probs * self._log(orig_probs))

        for kmer in all_kmers:
            information_gain = orig_entropy - self.conditional_entropy(kmer)
            if kwargs["threshold"] is not None and information_gain > kwargs["threshold"]:
                res[kmer] = information_gain
        
        return res


class KMerVariance(KMerInformationGain):
    """
    Select the k-mers that have high variance.
    """
    def __init__(self, kmer_set: KMerSet) -> None:
        super().__init__(kmer_set)
    
    def variance(self, kmer):
        positive_count = []
        negative_count = []

        for label in self.label_dict:
            # Count how many references contain this kmer
            positive_count.append(self.counter_dict[label][kmer])
            negative_count.append(self.label_dict[label] - self.counter_dict[label][kmer])
        
        positive_count = np.array(positive_count)
        negative_count = np.array(negative_count)

        p = positive_count / np.sum(positive_count + negative_count)
        return np.var(p)


    def select_kmers(self, **kwargs):
        res = dict()

        all_kmers = set()
        for label in self.label_dict:
            all_kmers.update(self.counter_dict[label].keys())

        for kmer in all_kmers:
            variance = self.variance(kmer)
            if kwargs["threshold"] is not None and variance > kwargs["threshold"]:
                res[kmer] = variance

        return res

class KMerChiSquaredTest(KMerInformationGain):
    """
    Do a Chi-squared test on the dependency of classification output on presence of k-mer.
    """
    def __init__(self, kmer_set: KMerSet) -> None:
        super().__init__(kmer_set)
    
    def chi_squared_test(self, kmer):
        from scipy.stats import chi2

        positive_count = []
        negative_count = []

        for label in self.label_dict:
            # Count how many references contain this kmer
            positive_count.append(self.counter_dict[label][kmer])
            negative_count.append(self.label_dict[label] - self.counter_dict[label][kmer])
        
        # Observed counts
        positive_count = np.array(positive_count)
        negative_count = np.array(negative_count)
        total_count = positive_count + negative_count

        # Expected Counts
        positive_rate = np.sum(positive_count) / np.sum(total_count)
        positive_expected = positive_rate * total_count
        negative_expected = (1-positive_rate) * total_count

        # Do chi-squared test
        degree_of_freedom = len(positive_count) - 1
        statistics = np.sum((positive_count - positive_expected) ** 2 / positive_expected + \
                     (negative_count - negative_expected) ** 2 / negative_expected)
        p_value = 1 - chi2(degree_of_freedom).cdf(statistics)

        return p_value


    def select_kmers(self, **kwargs):
        res = dict()

        all_kmers = set()
        for label in self.label_dict:
            all_kmers.update(self.counter_dict[label].keys())

        for kmer in all_kmers:
            p_value = self.chi_squared_test(kmer)
            if "threshold" not in kwargs or p_value <= kwargs["threshold"]:
                res[kmer] = p_value

        return res



if __name__ == "__main__":
    from db_sketching.kmer_set import FracMinHash

    