from Bio import SeqIO

class RandomReadDataset:
    def __init__(self, file_list, file_labels, read_length, random_errors=None, random_reverse_complement=False, debug=False):
        self.references = self._load_reference(file_list)
        

    def _load_reference(self, reference_list):
        res = []
        for file in reference_list:
            if file.endswith('.fna') or file.endswith('.fasta') or file.endswith('.fa'):
                file_type = "fasta"
            else:
                file_type = "fastq"
            for record in SeqIO.parse(file, file_type):
                res.append(record.seq)
        return res