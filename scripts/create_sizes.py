#!/usr/bin/env python3

"""
Created: February 10, 2022
Updated: August 23, 2023
Author(s): Todd Lenz, tlenz001@ucr.edu

Finds the length of all chromosomes in a FASTA file and outputs a
tab-delimited text file. The output will be stored in the same directory
as the input FASTA file.
"""


import sys
from Bio import SeqIO
import os


g = SeqIO.parse(open(sys.argv[1]),"fasta")
out = open("".join([os.path.dirname(os.path.abspath(sys.argv[1])), "/", 
                    os.path.splitext(os.path.basename(sys.argv[1]))[0], 
                    ".chrom.sizes"]),"w")
for i in g:
    out.write(i.id + "\t" + str(len(i)) + "\n")
out.close()
