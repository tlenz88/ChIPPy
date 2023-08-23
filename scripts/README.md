# ChIPPy scripts
There are several scripts used by ChIPPy to generate input files or intermediate files required by other tools within the pipeline.
##

## create_sizes.py

Finds the length of all chromosomes in a FASTA file and outputs a tab-delimited text file. The output will be stored in the same directory as the input FASTA file.

**Usage**:

```create_sizes.py genome.fasta```

## generate_chip_metadata.py

Generates the ChIP-seq metadata file required as input for the DiffBind R package. The output may be incorrect if replicates do not share a similar naming convention. This script is meant to be used as a component of the ChIPPy pipeline, but can be used independently if the sorted BAM file names contain the string '_sorted.bam' and the XLS peak files contain the string '_peaks.xls'.

**Usage**:
```generate_chip_metadata.py -i INPUT -t TREATMENT -c CONTROL```

- **Required arguments**:

   --input [-i]: Directory containing folders with ChIP-seq data. Each sample folder should contain a coordinate sorted BAM file and the XLS file with called peaks output by MACS.

   --treatment [-t]: Name or value indicating the targeted antibody.

   --control [-c]: Name or value indicating the contol antibody or input.

