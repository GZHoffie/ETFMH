import numpy as np

class MutationModel:
    def __init__(self) -> None:
        pass

    def containment(self, k):
        pass

class RandomMutationModel(MutationModel):
    """
    Assume mutation is distributed uniformly at random in the read, and
    all k-mers are unique. How will the containment index change?
    """
    def __init__(self, ANI, length) -> None:
        super().__init__()
        self.p = 1 - ANI
        self.ground_truth = np.random.choice(a=[False, True], size=length, p=[self.p, 1-self.p])
    
    def containment(self, k):
        kmer_matching = np.all(np.lib.stride_tricks.sliding_window_view(self.ground_truth, window_shape=k), axis=1)
        return kmer_matching.sum() / len(kmer_matching)


if __name__ == "__main__":
    m = RandomMutationModel(0.9, 1000)
    for i in range(1, 51):
        containment = m.containment(i)
        print(i, containment, containment ** (1/i))

        
