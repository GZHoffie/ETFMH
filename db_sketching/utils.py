from Bio.Seq import Seq
from tqdm import tqdm

class Seq2KMers:
    """
    Util class that outputs the sequence of K-mers given an input sequence.
    """

    def __init__(self, seed_length : int):
        """
        Args:
            - seed_length (int): length of k-mer.
        """
        # Initialize parameters for the minimizer.
        self._k = seed_length

        # Utils to find minimizers
        self._nucleotide_to_hash = {'A': 0, 'a': 0,
                                    'C': 1, 'c': 1,
                                    'G': 2, 'g': 2,
                                    'T': 3, 't': 3}
        self._seed_mask = int("1" * 2 * self._k, 2)


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
            
            hash_values.append(current_hash)
            
            # find the hash of rest of the k-mers
            for i in range(len(sequence) - self._k):
                current_hash = (current_hash << 2) & self._seed_mask
                current_hash = current_hash | sequence[i + self._k]
                hash_values.append(current_hash)
        
        return hash_values

    def canonical_kmers(self, sequence : str):
        """
        Given an input sequence, find all the canonical k-mers inside, in the form of a list of hash values.
        """
        forward_hash = self._hash(sequence)
        reverse_hash = self._rev_comp(forward_hash)

        forward_kmers = self._kmers(forward_hash)
        reverse_kmers = self._kmers(reverse_hash)

        min_kmers = [min(i, j) for i, j in zip(forward_kmers, reversed(reverse_kmers))]
        return min_kmers



if __name__ == "__main__":
    seq2kmers = Seq2KMers(3)
    print(seq2kmers.canonical_kmers("ACGTGGGCGT"))


