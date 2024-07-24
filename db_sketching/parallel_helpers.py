import sys
import subprocess
from Bio.SeqIO import read, parse
import orthoani
from multiprocess import Pool

sys.path.append('../')
sys.path.append('../../')

from db_sketching.kmer_set import EstimatorFracMinHash

    
def check_genome_files(genome_files):
    checked_genome_files = []
    for file in genome_files:
        try:
            parsed_file = parse(file,"fasta")
            assert(len([record for record in parsed_file]) > 0)
            checked_genome_files.append(file)
        except:
            print(f"File {file} is damaged / invalid")

    return checked_genome_files

def get_genome_length(genome_file):
    parsed_file = parse(genome_file,"fasta")
    return sum(len(record) for record in parsed_file)


def compute_length_parallel(genome_files,max_pool=16):
    with Pool(max_pool) as p:
        return p.starmap(get_genome_length,([(g_file,) for g_file in genome_files]))
    
def compute_ortho_ani(genome_file_1, genome_file_2):
    try:
        genome_1_read = parse(genome_file_1,"fasta")
        genome_2_read = parse(genome_file_2,"fasta")
        ortho_ani_value = orthoani.orthoani(genome_1_read,genome_2_read)
        return ortho_ani_value
    except:
        return 0

def compute_ortho_ani_parallel(genome_files_1,genome_files_2,max_pool=16):
    args = [(g1,g2) for g1,g2 in zip(genome_files_1,genome_files_2)]
    with Pool(max_pool) as p:
        return p.starmap(compute_ortho_ani,args)

def compute_pairwise_ortho(genome_files):
    genome_files_1, genome_files_2 = zip(*[(g1,g2) for g1 in genome_files for g2 in genome_files])
    return compute_ortho_ani_parallel(genome_files_1,genome_files_2)

def compute_estimator_kmer_sketches(genome_file, conditions, kmer_template, canonical, multiplicity, estimator):
    genome_kmer = EstimatorFracMinHash(
        conditions=conditions, 
        kmer_template=kmer_template, 
        canonical=canonical, 
        multiplicity=multiplicity, 
        estimator=estimator
    )
    genome_kmer.insert_file(genome_file)
    return genome_kmer

def compute_estimator_kmer_sketches_parallel(genome_files, conditions, kmer_template, canonical, multiplicity, estimator,max_pool=16):
    args = [(g,conditions, kmer_template, canonical, multiplicity, estimator) for g in genome_files]
    with Pool(max_pool) as p:
        return p.starmap(compute_estimator_kmer_sketches,args)


def compute_kmer_ani(genome_1_kmer, genome_2_kmer):    
    kmer_estimated_ani = genome_1_kmer.ANI_estimation(genome_2_kmer)
    return kmer_estimated_ani


def compute_kmer_ani_parallel(kmer_sketches_1,kmer_sketches_2,max_pool=16):
    args = [(s1,s2) for s1,s2 in zip(kmer_sketches_1,kmer_sketches_2)]
    with Pool(max_pool) as p:
        return p.starmap(compute_kmer_ani,args)