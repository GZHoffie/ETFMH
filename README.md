# Experiments on ANI Estimation

The main files are under `db_sketching` directory. The `utils.py` file contains a utility function that converts a DNA sequence into a list of k-mer hash values, and `kmer_set.py` is a basic implementation of k-mer set, as well as functions to calculate resemblence, containment index, and estimation of ANI.

## Installation

To use the experimental code, you need to insall biopython,

```shell
pip install biopython
```

and add the python package into `PYTHONPATH`,

```shell
git clone https://github.com/GZHoffie/ETFMH.git
cd ETFMH
export PYTHONPATH=$(pwd)
```

## Usage

Then, we can use it to estimate ANI between two reference genomes in `./data/`.

```python
from db_sketching.kmer_set import KMerSet, FracMinHash
from Bio import SeqIO

# Path to the datasets we are testing
data_1 = ["../data/562.fna"]
data_2 = ["../data/564.fna"]

# The boolean function used in FracMinHash
def cond(kmer_hash):
    hash = (976369 * kmer_hash + 1982627) % 10000
    if hash < 100:
        return True
    else:
        return False

   
frac1 = FracMinHash(cond, k)
frac2 = FracMinHash(cond, k)

# insert dataset
frac1.insert_file_list(data_1)
frac2.insert_file_list(data_2)

# Calculate containment index and estimate ANI
print(frac1.containment(frac2))
print(frac1.ANI_estimation(frac2))
```

Refer to the python notebook in `./test` for details.
