#!/usr/bin/env python3

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
