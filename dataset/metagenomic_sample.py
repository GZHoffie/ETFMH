import pandas as pd
import numpy as np
import subprocess
import os
import random
import glob
from pathlib import Path

class MetagenomicSampleGenerator:
    def __init__(self) -> None:
        pass

    def generate_sample(self, references, read_num, output_directory, output_file_name, distribution="log_normal"):
        """
        Generate a metagenomic sample.

        Install badread with
        ```bash
        pip3 install git+https://github.com/rrwick/Badread.git
        ```
        before running this function.

        Args:
            - read_num (int): number of reads in the generated sample.
            - references (List[str]): list of file names containing the references in the metagenomic sample.
            - output_directory (str): path to store the output file.
            - output_file_name (str): name of the output files.
            - distribution (either "uniform" or "log_normal"): the distribution of abundance.
        """
        # decide abundance of the sample
        total_species_num = len(references)
        if distribution == "log_normal":
            species_abundance = np.random.lognormal(0, 1, total_species_num)
            species_abundance = species_abundance / np.sum(species_abundance)
            species_abundance = np.sort(species_abundance)
        else:
            species_abundance = np.ones(total_species_num) / len(total_species_num)


        # Assign pathogens the lowest abundance
        species_index = 0
        reads = []
        all_species = []
        abundance_ground_truth = {}

        # Simulate reads from pathogen genome
        #print(references)
        for file in references:
            # determine number of reads of this species
            num_reads = int(species_abundance[species_index] * read_num * 4000)

            # Simulate reads
            simulated_read = subprocess.run(["badread", "simulate", "--reference", file, "--quantity", str(num_reads), "--length", "4000,2000", 
                                             "--glitches", "0,0,0", "--junk_reads", "0", "--random_reads", "0", "--chimeras", "0"], capture_output=True)
            
            # Store the reads in `reads` list
            read_list = simulated_read.stdout.decode("utf-8").split('\n')
            read_list = [r for i, r in enumerate(read_list) if i % 4 == 1]
            reads.extend(read_list)

            # Accession = file.stem()
            accession = Path(file).stem

            # store ground truth
            abundance_ground_truth[accession] = species_abundance[species_index]

            # move on to next species
            species_index += 1

        # Output raw reads
        read_index = np.arange(len(reads))
        random.shuffle(read_index)

        read_i = 0
        with open(os.path.join(output_directory, output_file_name + ".fasta"), 'w') as f:
            for i, r in enumerate(read_index):
                f.write(">" + str(i) + "\n")
                f.write(reads[r] + "\n")


        # Output ground truth abundance
        with open(os.path.join(output_directory, output_file_name + "_abundance.csv"), 'w') as f:
            f.write("accession,abundance\n")
            for file in abundance_ground_truth:
                f.write(file + "," + str(abundance_ground_truth[file]) + "\n")


if __name__ == "__main__":
    import glob

    g = MetagenomicSampleGenerator()
    #print(glob.glob("./data_temp/*/*.fna"))
    g.generate_sample(glob.glob("./single_species/*/*.fna"), 100000, "./", "EColi")