import glob
import pathlib
import numpy as np
import shutil
import os
import subprocess

def split_dataset(directory_name, target_directory, training_set_proportion=0.8):
    items = glob.glob(directory_name + "*")
    training_items = np.random.choice(items, int(training_set_proportion * len(items)), replace=False)
    test_items = [i for i in items if i not in training_items]

    # Put the training and test items into respective directory
    subprocess.run(["mkdir", "-p", target_directory + "/train"])
    subprocess.run(["mkdir", "-p", target_directory + "/test"])
    for i in training_items:
        subprocess.run(["cp", "-r", i, target_directory + "/train/"])
    
    for i in test_items:
        subprocess.run(["cp", "-r", i, target_directory + "/test/"])



if __name__ == "__main__":
    split_dataset("/home/zhenhao/ETFMH/data/Staphylococcus/", "/home/zhenhao/ETFMH/data/Staphylococcus_data/")
