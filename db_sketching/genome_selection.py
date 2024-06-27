from db_sketching.kmer_set import KMerSet
import numpy as np
import copy
from unionfind import unionfind

class GenomeFiltering:
    def __init__(self, kmer_set : KMerSet) -> None:
        pass
        self.genome_dict = dict()
        self.kmer_set = kmer_set
    
    def insert_genome(self, genome_file_name):
        # insert the new genome
        self.kmer_set.insert_file(genome_file_name)

        # Store the kmer set in self.genome_dict
        self.genome_dict[genome_file_name] = copy.deepcopy(self.kmer_set)
        self.kmer_set.reset()
    
    def _calculate_pairwise_distance(self):
        """
        Calculate the pairwise distance using Jaccard index, and store it in a 
        distance matrix.
        """
        genome_list = list(self.genome_dict.keys())
        distance_matrix = np.zeros((len(genome_list), len(genome_list)))
        for i in range(len(genome_list)):
            for j in range(i+1, len(genome_list)):
                distance = self.genome_dict[genome_list[i]].resemblence(self.genome_dict[genome_list[j]])
                distance_matrix[i][j] = distance_matrix[j][i] = distance
        
        return genome_list, distance_matrix
    
    def hierarchical_clustering(self, genome_list, distance_matrix, threshold):
        u = unionfind(len(genome_list))
        for i in range(len(genome_list)):
            for j in range(i+1, len(genome_list)):
                if distance_matrix[i][j] <= threshold:
                    u.unite(i, j)
        
        return u.groups()


