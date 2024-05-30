from Bio.Seq import Seq
from tqdm import tqdm
import heapdict

class Seq2KMers:
    """
    An implementation of the Minimizer algorithm

    Args:
        vocab (:obj:`int`, `optional`):
            A dictionnary of k-mer keys and their ids :obj:`{AAA (0): 0,...}`

        unk_token (:obj:`str`, `optional`):
            The unknown token to be used by the model.
    """

    def __init__(self, seed_length, window_length=None, use_mask=True):
        # Initialize parameters for the minimizer.
        self._k = seed_length
        self._l = window_length if window_length is not None else seed_length

        if (self._l < self._k):
            print(f"[ERROR]\t\tthe window length (currently {self._l}) must be equal or larger than the seed length (currently {self._k}).")
            return

        # Utils to find minimizers
        self._nucleotide_to_hash = {'A': 0, 'a': 0,
                                    'C': 1, 'c': 1,
                                    'G': 2, 'g': 2,
                                    'T': 3, 't': 3}
        self._seed_mask = int("1" * 2 * self._k, 2)
        # Random permuatation of hash values, same as SeqAn3
        self._hash_mask = 0x8F3F73B5CF1C9ADE if use_mask else 0


    def _hash(self, sequence):
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

    def canonical_kmers(self, sequence):
        forward_hash = self._hash(sequence)
        reverse_hash = self._rev_comp(forward_hash)

        forward_kmers = self._kmers(forward_hash)
        reverse_kmers = self._kmers(reverse_hash)

        min_kmers = [min(i, j) for i, j in zip(forward_kmers, reversed(reverse_kmers))]
        return min_kmers


    def minimizers(self, sequence):
        """
        Find the minimizers of the sequence.
        """
        minimizers = []
        if len(sequence) - self._k + 1 <= 0:
            return minimizers
        
        # Find hash values of k-mers
        forward_hash = self._hash(sequence)
        reverse_hash = self._rev_comp(forward_hash)

        masked_forward_kmers = [h ^ self._hash_mask for h in self._kmers(forward_hash)]
        masked_reverse_kmers = [h ^ self._hash_mask for h in self._kmers(reverse_hash)]

        # find canonical k-mers
        masked_min_kmers = [min(i, j) for i, j in zip(masked_forward_kmers, reversed(masked_reverse_kmers))]

        # Find min in sliding window
        sliding_window = heapdict.heapdict()
        for i in range(min(self._l - self._k + 1, len(masked_min_kmers))):
            sliding_window[i] = masked_min_kmers[i]

        if len(sliding_window) != 0:
            minimizers.append(sliding_window.peekitem()[1])

        for i in range(len(masked_min_kmers) - self._l):
            sliding_window.pop(i)
            current_index = i + self._l - self._k + 1
            sliding_window[current_index] = masked_min_kmers[current_index]
            if len(sliding_window) != 0:
                minimizer = sliding_window.peekitem()[1]
                if minimizers[-1] != minimizer:
                    minimizers.append(minimizer)

        unmasked_minimizers = [m ^ self._hash_mask for m in minimizers]
        tokens = [str(m) for m in unmasked_minimizers]
        #print(tokens)

        return tokens
    
    def init_vocab(self, condition):
        MASK = 3 # Mask to filter out single nucleotide in k-mer
        vocab = {}
        vocab_index = 0
        for i in tqdm(range(4 ** self._k), desc=f"[INFO]\t\tBuilding up vocabulary using {self._k}-mers"):
            kmer_sequence = [(i >> (j * 2)) & MASK for j in reversed(range(self._k))]
            rev_comp_sequence = self._rev_comp(kmer_sequence)

            kmer_hash = 0
            rev_comp_hash = 0
            for i in range(self._k):
                kmer_hash = kmer_hash << 2
                rev_comp_hash = rev_comp_hash << 2
                kmer_hash = kmer_hash | kmer_sequence[i]
                rev_comp_hash = rev_comp_hash | rev_comp_sequence[i]
        
            # Pick the canonical k-mer only
            min_hash = kmer_hash if kmer_hash < rev_comp_hash else rev_comp_hash
            if condition(min_hash) and min_hash not in vocab:
                vocab[min_hash] = vocab_index
                vocab_index += 1
        

        #print(vocab)
        print(f"[INFO]\t\tVocabulary build complete. Vocab size: {len(vocab)}.")
        return vocab, len(vocab)


