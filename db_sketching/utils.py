from Bio.Seq import Seq
from tqdm import tqdm
import heapdict
import copy

class Seq2KMers:
    """
    Util class that outputs the sequence of K-mers given an input sequence.
    """

    def __init__(self, seed_template):
        """
        Args:
            - seed_template (str): a string consisting of 0 and 1's, indicating which base to include in the seed.
              For example, seed_template = "1" * 12 represents a contiguous 12-mer, "1001001" means taking one for every
              3 bases. 
        """
        # Initialize parameters for the minimizer.
        self._k = len(seed_template)

        # Utils to find minimizers
        self._nucleotide_to_hash = {'A': 0, 'a': 0,
                                    'C': 1, 'c': 1,
                                    'G': 2, 'g': 2,
                                    'T': 3, 't': 3}
        self._seed_mask = int("1" * 2 * self._k, 2)
        self._seed_template = int(''.join([char * 2 for char in seed_template]), 2)


    def _hash(self, sequence : str):
        """
        Convert a DNA sequence to a sequence of numbers in 0-3.
        """
        res = []
        temp = [] # storing everything before an ambiguous nucleotide
        for n in sequence:
            if n in self._nucleotide_to_hash:
                temp.append(self._nucleotide_to_hash[n])
            elif len(temp) > 0:
                res.append(temp)
                temp = []
        
        if len(temp) > 0:
            res.append(temp)
        
        return res
    

    def _rev_comp(self, sequence_hash):
        """
        Given a sequence of hash values returned from self._hash(sequence),
        return the reverse complement of that sequence.
        """
        res = []
        for sequence in reversed(sequence_hash):
            res.append([(3 - h) for h in reversed(sequence)])

        return res


    def _kmers(self, sequence_hash):
        """
        Convert k-mers to their hash values.
        """
        hash_values = []
        for sequence in sequence_hash:
            if len(sequence) - self._k + 1 <= 0:
                continue
            
            # find the hash value of the first k-mer
            current_hash = 0
            for i in range(self._k):
                current_hash = current_hash << 2
                current_hash = current_hash | sequence[i]
            
            hash_values.append(current_hash & self._seed_template)
            
            # find the hash of rest of the k-mers
            for i in range(len(sequence) - self._k):
                current_hash = (current_hash << 2) & self._seed_mask
                current_hash = current_hash | sequence[i + self._k]
                hash_values.append(current_hash & self._seed_template)
        
        return hash_values
    
    def kmers(self, sequence : str):
        return self._kmers(self._hash(sequence))

    def canonical_kmers(self, sequence : str):
        """
        Given an input sequence, find all the canonical k-mers inside, in the form of a list of hash values.
        NOTE: canonical k-mers with spaced seeds is not well defined. Use this only when using contiguous seeds.
        """
        forward_hash = self._hash(sequence)
        reverse_hash = self._rev_comp(forward_hash)

        forward_kmers = self._kmers(forward_hash)
        reverse_kmers = self._kmers(reverse_hash)

        min_kmers = [min(i, j) for i, j in zip(forward_kmers, reversed(reverse_kmers))]
        return min_kmers

    def minimizers(self, sequence, window_size, hash_mask=0):
        """
        Find the minimizers of the sequence.
        NOTE: canonical k-mers with spaced seeds is not well defined. Use this only when using contiguous seeds.

        Args:
            - sequence (str): the sequence to be converted.
            - window_size (int): The window size in which to find minimizer.
            - hash_mask (int): Can be a large integer to indicate a random permutation
              of k-mer hash values.
        """
        minimizers = []
        if len(sequence) - self._k + 1 <= 0:
            return minimizers
        
        if window_size < self._k:
            print("[ERROR]\t\tThe window size must be at least k.")
            return minimizers
        
        # Find hash values of k-mers
        forward_hash = self._hash(sequence)
        reverse_hash = self._rev_comp(forward_hash)

        masked_forward_kmers = [h ^ hash_mask for h in self._kmers(forward_hash)]
        masked_reverse_kmers = [h ^ hash_mask for h in self._kmers(reverse_hash)]

        # find canonical k-mers
        masked_min_kmers = [min(i, j) for i, j in zip(masked_forward_kmers, reversed(masked_reverse_kmers))]

        # Find min in sliding window
        sliding_window = heapdict.heapdict()
        for i in range(min(window_size - self._k + 1, len(masked_min_kmers))):
            sliding_window[i] = masked_min_kmers[i]

        if len(sliding_window) != 0:
            minimizers.append(sliding_window.peekitem()[1])

        for i in range(len(masked_min_kmers) - window_size):
            sliding_window.pop(i)
            current_index = i + window_size - self._k + 1
            sliding_window[current_index] = masked_min_kmers[current_index]
            if len(sliding_window) != 0:
                minimizer = sliding_window.peekitem()[1]
                if minimizers[-1] != minimizer:
                    minimizers.append(minimizer)

        unmasked_minimizers = [m ^ hash_mask for m in minimizers]
        tokens = [m for m in unmasked_minimizers]
        #print(tokens)

        return tokens
    
    def all_canonical_kmers(self):
        """
        Return a dictionary that maps all k-mers to their canonical k-mers.
        NOTE: canonical k-mers with spaced seeds is not well defined. Use this only when using contiguous seeds.
        """
        MASK = 3 # Mask to filter out single nucleotide in k-mer
        vocab = {}
        vocab_index = 0
        for i in tqdm(range(4 ** self._k), desc=f"[INFO]\t\tBuilding up vocabulary using {self._k}-mers"):
            if i in vocab:
                continue

            kmer_sequence = [[(i >> (j * 2)) & MASK for j in reversed(range(self._k))]]
            rev_comp_sequence = self._rev_comp(kmer_sequence)

            kmer_hash = 0
            rev_comp_hash = 0
            for i in range(self._k):
                kmer_hash = kmer_hash << 2
                rev_comp_hash = rev_comp_hash << 2
                kmer_hash = kmer_hash | kmer_sequence[0][i]
                rev_comp_hash = rev_comp_hash | rev_comp_sequence[0][i]
        
            # Pick the canonical k-mer only
            min_hash = kmer_hash if kmer_hash < rev_comp_hash else rev_comp_hash
            vocab[kmer_hash] = min_hash
            vocab[rev_comp_hash] = min_hash
        

        #print(vocab)
        print(f"[INFO]\t\tNumber of canonical k-mers: {len(set(vocab.values()))}.")
        return vocab


class KMer:
    """
    A util class to find all 
    """
    def __init__(self, k) -> None:
        self._k = k
        self.kmer = None
    
    def kmer_to_list(self, kmer: int):
        """
        Represent the given kmer in list form.
        """
        return [(kmer >> (j * 2)) & 3 for j in reversed(range(self._k))]
    
    def list_to_kmer(self, kmer_list):
        """
        Given kmer in the form of a list of integers, transform it back to integer form.
        """
        res = 0
        for i in range(self._k):
            res = res << 2
            res = res | kmer_list[i]
        return res

    
    def reverse_complement(self, kmer_list):
        rev_comp_sequence = [(3 - h) for h in reversed(kmer_list)]
        #rev_comp_hash = self.list_to_kmer(rev_comp_sequence)
        
        return rev_comp_sequence
    
    def canonical_kmer(self, kmer: int):
        """
        Given a k-mer, return the corresponding canonical k-mer.
        """
        rev_comp_sequence = self.reverse_complement(self.kmer_to_list(kmer))
        rev_comp_hash = self.list_to_kmer(rev_comp_sequence)
        return min(rev_comp_hash, kmer)

    
    def distance_one_neighbors(self, kmer: int, distance: str = "hamming", canonical: bool = True):
        """
        Return a list of hash values that has distance one to self.kmer.

        Args:
            - distance (str): either "hamming" or "edit", indicating hamming distance or edit distance.
            - canonical (bool): whether we include the reverse complement of the neighbors as well.
        """
        neighbors_list = list()

        # Represent the original kmer in list form
        orig_kmer = self.kmer_to_list(kmer)

        # Include all kmers with hamming distance 1
        for i in range(self._k):
            for j in range(4):
                modified_kmer = copy.deepcopy(orig_kmer)
                modified_kmer[i] = j
                neighbors_list.append(modified_kmer)
        
        # TODO: implement the case for edit distance
        if canonical:
            neighbors_list_rev_comp = list()
            for i in neighbors_list:
                neighbors_list_rev_comp.append(self.reverse_complement(i))
            
            neighbors_list += neighbors_list_rev_comp
        
        
        return set([self.list_to_kmer(i) for i in neighbors_list])



if __name__ == "__main__":
    seq2kmers = Seq2KMers("101")
    #print(seq2kmers.all_canonical_kmers())
    print(seq2kmers.kmers("ACGTGGGCGT"))

    #a = KMer(19)
    #print(a.distance_one_neighbors(0))


