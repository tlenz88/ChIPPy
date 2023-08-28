# ChIPPy Scripts
There are several scripts used by ChIPPy to generate input files or intermediate files required by other tools within the pipeline.

## create_sizes.py

Finds the length of all chromosomes in a FASTA file and outputs a tab-delimited text file. The output will be stored in the same directory as the input FASTA file.

**Usage**: ```python3 create_sizes.py *.fasta```

## differential_peak_calling.R

Performs differential peak calling using a metadata file describing the experimental setup. See ```example_metadata.txt``` in the ```examples``` folder. 

**Usage**: ```Rscript differential_peak_calling.R chip_metadata.txt```

## peaks2genes.py

Finds peaks identified by ```macs2 callpeak``` whose summit is contained within genes and 5' promoter regions in a ```gff``` file. Accepts either ```narrowPeak``` or ```broadPeak``` input.

**Usage**: ```python3 peaks2genes.py *.gff *.narrowPeak```
