#!/usr/bin/env python
from Bio import SeqIO
import math
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", dest = "infile", help = "The fasta-file which contigs you want to rename")
parser.add_argument("-o", "--output_dest", dest = "outfile", help = "The file path where you want the output (will overwrite things)")
args = parser.parse_args()
binfile = args.infile
boutfile = args.outfile

with open(binfile, "rt") as bin:
    all_contigs = list(SeqIO.parse(bin, "fasta"))
    num_len = int(math.log10(len(all_contigs))+1)
    bin_name = binfile.split('.',1)[0]
    for i,s in enumerate(all_contigs):
        s.id =f"{bin_name}_{i+1:0{num_len}}" #renames contigs
    with open(boutfile, "w") as output_handle:
        SeqIO.write(all_contigs, output_handle, "fasta")
