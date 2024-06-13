import numpy as np
import matplotlib.pyplot as plt
import statistics
import random

class MutationModel:
    """
    Base class for a mutation model 

    Contains a ground truth 
    """
    def __init__(self) -> None:
        pass

    def containment(self, k):
        pass

class TrueFalseMutationModel(MutationModel):
    """
    Ground truth represented as a numpy array of True/False values 
    Containment computed by taking k-mers across the ground truth using a sliding window
    """
    def containment(self, k):
        kmer_matching = np.all(np.lib.stride_tricks.sliding_window_view(self.ground_truth, window_shape=k), axis=1)
        return kmer_matching.sum() / len(kmer_matching)

class RandomMutationModel(TrueFalseMutationModel):
    """
    Assume mutation is distributed uniformly at random in the read, and
    all k-mers are unique. How will the containment index change?
    """
    def __init__(self, ANI, length) -> None:
        super().__init__()
        self.p = 1 - ANI
        self.ground_truth = np.random.choice(a=[False, True], size=length, p=[self.p, 1-self.p])
    

class Prev1MutationModel(TrueFalseMutationModel):
    """
    p0 : Probability that the current base will mutate if the previous base does not mutate
    p1 : Probability that the current base will not mutate if the previous base mutated

    In steady state, ANI = (1 - p1) / (1 - p1 + p0)
    """
    def __init__(self, p0, p1, length) -> None:
        super().__init__()
        self.ani = (1 - p1) / (1 - p1 + p0)

        prev_mut = True
        self.ground_truth = []
        for _ in range(length):
            if prev_mut:
                self.ground_truth.append(np.random.choice(a=[False, True],p=[p0, 1-p0]))
            else:
                self.ground_truth.append(np.random.choice(a=[False, True],p=[p1, 1-p1]))
            prev_mut = self.ground_truth[-1]

        self.ground_truth = np.array(self.ground_truth)
        self.ground_ani = (self.ground_truth * 1).sum() / len(self.ground_truth)



if __name__ == "__main__":
    ANI = 0.90
    LENGTH = 10000
    KMAX = 40
    SAMPLES = 1000

    # m = RandomMutationModel(ANI, LENGTH)
    # for i in range(1, KMAX + 1):
    #     containment = m.containment(i)
    #     print(i, containment, containment ** (1/i))

    