import glob
import pathlib
import numpy as np
import shutil
import os
import subprocess

def split_dataset(directory_name, target_directory, training_set_species_proportion=0.8, training_set_strains_proportion=0.8):
    species = glob.glob(directory_name + "*")
    training_species = np.random.choice(species, int(training_set_species_proportion * len(species)), replace=False)

    test_species = [i for i in species if i not in training_species]

    # Put the training and test items into respective directory
    subprocess.run(["mkdir", "-p", target_directory + "/train"])
    subprocess.run(["mkdir", "-p", target_directory + "/test_ood_species"])
    subprocess.run(["mkdir", "-p", target_directory + "/test_ood_strains"])
    for i in training_species:
        subprocess.run(["cp", "-r", i, target_directory + "/train/"])
    
    # out-of-domain strains
    for i in glob.glob(target_directory + "/train/*"):
        species_name = i.split("/")[-1]
        subprocess.run(["mkdir", "-p", target_directory + "/test_ood_strains/" + species_name])
        strains = glob.glob(i + "/*.fna")
        test_strains = np.random.choice(strains, int((1-training_set_strains_proportion) * len(strains)), replace=False)
        for j in test_strains:
            subprocess.run(["mv", j, target_directory + "/test_ood_strains/" + species_name])

    # out-of-domain species
    for i in test_species:
        subprocess.run(["cp", "-r", i, target_directory + "/test_ood_species/"])



if __name__ == "__main__":
    split_dataset("/home/zhenhao/ETFMH/data/Staphylococcus/", "/home/zhenhao/ETFMH/data/Staphylococcus_data/")
