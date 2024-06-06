from pathlib import Path
import subprocess
import os
import glob
from tqdm import tqdm
import pandas as pd
import pickle
from multiprocessing import Process, Queue, Pool



class ReferenceGenomeDownloader:
    """
    A utility class enabling us to download reference genomes from NCBI RefSeq database.
    """
    def __init__(self, metadata_df) -> None:
        """
        Args:
            - metadata_df (pd.DataFrame): A dataframe consisting of GTDB metadata.
        """
        self.metadata_df = metadata_df
    
    def _download_reference(self, accession, directory, name):
        """
        Download the reference using the given accession number, and store it in
        ${directory}/${name}.fna. If the file already exists, we append to the file.
        """
        # mkdir if the directory doesn't exist
        Path(directory).mkdir(parents=True, exist_ok=True)

        # Change to this directory
        os.chdir(directory)

        # Download using NCBI datasets
        subprocess.run(["datasets", "download", "genome", "accession",
                         accession, "--filename", name + ".zip"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if os.path.exists(name + ".zip"):
            subprocess.run(["unzip", name + ".zip"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Store the results in ${directory}/${taxid}.fna
        files_name = glob.glob("ncbi_dataset/data/*/*.fna")
        with open(name + ".fna", "wb") as f:
            subprocess.run(["cat"] + files_name, stdout=f)
        subprocess.run(["rm", "-r", "ncbi_dataset", name + ".zip", "README.md"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def _remove_directory(self, directory):
        """
        Delete the directory and all the content inside.
        """
        subprocess.run(["rm", "-r", directory], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


    def find_signature(self, accession_list, taxid_list, tax_id):
        """
        Given a list of accession, determine the signature.
        """
        # Download the reference corresponding to the accession
        directory = "/home/zhenhao/TDT/data_temp/" + str(tax_id)

        for i, accession in tqdm(enumerate(accession_list), desc="Downloading references for genus " + str(tax_id)):
            self._download_reference(accession, directory, str(i))

        #print(glob.glob(directory + "/*.fna"))

        res = self.signature.find_consensus(files=glob.glob(directory + "/*.fna"))
        self._remove_directory(directory)

        return tax_id, res
    
    def find_all_signatures(self, metadata_df, num_samples=None):
        """
        Given a metadata_df with columns "species_taxid", "genus_taxid", "ncbi_genbank_assembly_accession"
        find the signature for all genera.

        For each genera, only take `num_samples` species to speed things up.
        """
        res = {}
        family_set = set(metadata_df.dropna(subset=['family_taxid'])['family_taxid'])
        argument_list = []
        for i, family in enumerate(family_set):
            print("Processing family", i+1, "/", len(family_set))
            # sample just one genome per species
            family_df = metadata_df[metadata_df["family_taxid"] == family].groupby("species_name").sample(1)

            # if there are more than `num_samples` species, sample only `num_samples` of them
            if num_samples is not None and len(family_df) > num_samples:
                family_df = family_df.sample(num_samples)
            
            # Download the references
            argument_list.append((family_df["ncbi_genbank_assembly_accession"],
                                  family_df["species_taxid"], family))
            #signature = self.find_signature(genus_df["ncbi_genbank_assembly_accession"],
            #                                genus_df["species_taxid"], genus)
            #res[genus] = signature
        #return res

        # Use multiprocessing
        print("Using", os.cpu_count(), "CPUs.")
        pool = Pool(os.cpu_count())
        res = pool.starmap(self.find_signature, argument_list)
        
        return res
