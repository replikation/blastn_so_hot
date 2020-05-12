#!/usr/bin/env python
from Bio.Blast import NCBIWWW
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='ncbi based blast')

parser.add_argument(
    '--fasta', required=True,
    help='fasta file')

args = parser.parse_args()
fasta_file = args.fasta


record = SeqIO.read(fasta_file, format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"), results_file='my_blast.xml')

with open("my_blast.xml", "w") as out_handle:
     out_handle.write(result_handle.read())
result_handle.close()

