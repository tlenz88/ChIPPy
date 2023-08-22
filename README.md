# ChIPpipe
##
This is a complete pipeline for ChIP-seq data analysis. Most of the required packages are installed automatically. Pipeline can be started from any step, but specific files are required for the desired step.

## How do I organize my data?

Put all input files in a single directory with one folder per sample. The sequences can be either single or paired-end and can be gzipped, though this is not a requirement. Paired-end files need to have '_R1' and '_R2' within the file name to specify forward and reverse reads.

```
    + PATH_TO_DATA
        + sample1
            ++ sample1_R1.fastq.gz
            ++ sample1_R2.fastq.gz
        + sample2
            ++ sample2_R1.fastq.gz
            ++ sample2_R2.fastq.gz
        + sample 3
            ++ sample3.fastq.gz
        + sample 4
            ++ sample4.fastq.gz
```

## What tools do I need to run the pipeline?

Although the ChIPpipe script automatically downloads most of the necessary software/tools, depending on your system you may need to install some things manually. The only software that requires manual installation is [**Trimmomatic**](http://www.usadellab.org/cms/?page=trimmomatic) due to compatibility issues with the java version specified in the Trimmomatic build.xml file. Here is a list of all the necessary tools to run the complete pipeline:

- Python3/pip3
- R
- fastqc
- Trimmomatic
- Bowtie2
- Picardtools
- Samtools
- MACS2

Any version of Python3 and R are acceptable, but it is recommended that you use upgrade to the latest versions.

## What else do I need?

Any genome can be used as long as a properly formatted '.fasta' file is provided. For the tools that require a '.gff' file, ensure that the basename of the file is the same as the '.fasta' file provided. Any additional files, such as bowtie2 indexes and a chromosome sizes file, will be created automatically if not already present.

You may want to download IGV because a '.wig' file will be generated for easier visualization in the IGV genome browser.

## How do I use it?

Brief description of input arguments via help message:

```
    ChIPpipe.sh --help
    usage : ChIPpipe.sh -i INPUT -g GENOME [-o OUTPUT] [-s STEP] [-q QUALITY] [-t TREATMENT] [-c CONTROL] [-p THREADS] [-r]

    ----------------------------------------------------------------
    Required inputs:
      -i|--input  INPUT        : Input data folder.
      -g|--genome GENOME       : Path to genome files.

    Optional inputs:
      -o|--output OUTPUT       : Output folder.
      -s|--step STEP           : Choose starting step.
            quality_check      : Initial quality check.
                 trimming      : Adapter trimming.
                alignment      : Read alignment.
            deduplication      : Remove PCR duplicates.
                filtering      : Filtering low-quality reads.
                  sorting      : Sorting reads by coordinate.
                  mapping      : Mapping reads to each base pair
      -q|--quality QUALITY     : Phred quality score for filtering.
      -t|--treatment TREATMENT : Target antibody.
      -c|--control CONTROL     : Control antibody (IGG) or input.
      -p|--proc THREADS        : Processor threads.
      -r|--remove REMOVE       : Remove intermediate files.
      -h|--help HELP           : Show help message.
    ----------------------------------------------------------------
```

The following is a more detailed description of input arguments:

- **Required arguments**:

    --input [-i]: Directory containing a folder for each sample. Samples can be single or paired-end, but paired-end sequences need to contain '_R1' and '_R2' within the file names to determine forward and reverse reads. If starting from a step in the pipeline after alignment the input files are not required to contain '_R1' or '_R2'.

    --genome [-g]: Directory containing genome files for the organism of interest. A '.fasta' file should be present in the directory so that any additional required genome files can be automatically created if missing. All genome files should have the same basename.

- **Optional arguments**:

    --output [-o]: Output directory. If no output is given, the output files will be saved to the input folders.

    --step [-s]: Desired step at which to start the pipeline. If no step is entered, the entire pipeline will be run.

    --quality [-q]: Integer indicating Phred quality score for trimming low-quality bases at the ends of reads during read pairing and for removing low-quality reads during filtering. Phred score will also be used to approximate filtering settings during alignment.

    --treatment [-t]: Name or value indicating the targeted antibody in your experiment. This value is used to differentiate the targeted/treatment samples from the control samples, therefore the argument should be a string in the filenames. e.g. ```-t H3K9me3``` would indicate that all samples with 'H3K9me3' as part of their filename are the targeted/treatment samples. This arguement does not have to be valid antibody, but can be any string that is used within your sample names to indicate your targeted samples. However, if you are targeting a histone modification, this arguement should at least be the base name of your targeted histone modification, such as ```-t H3K``` or ```-t H2A```, so that peak calling is made more accurate by using the proper peak calling algorithm.

    --control [-c]: Name or value indicating the control or input samples in your experiment. Like '--antibody1', this value is used to differentiate the control samples from the targeted samples. This can also be any string that is used within your sample names, but will be used to indicate your control samples. e.g. ```-a2 IGG``` would indicate that all samples with 'IGG' as part of their filename are control samples.

    --proc [-p]: Integer indicating number of processor threads to use for tools that allow multithreading. If no value is given, the number of available threads will be determined automatically and will use half of available threads on your system.

    --remove [-r]: Removes intermediate files when no longer needed by pipeline to save hard drive space.

## What are the individual steps within the pipeline?

- **Quality check**:

    FastQC is used to check sequence quality prior to running analysis pipeline. In brief, statistics such as per base sequence quality, per base sequence content, sequence duplication levels and adapter content are calculated from your input '.fastq(.gz)' files. See [project website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for more details.

- **Trimming**:

    Trimmomatic is used to perform initial trimming. Adapter sequences in the provided 'adapters.txt' file and bases below the quality threshold set by the 'quality [-q]' argument (default = 30) are removed from the ends of reads. The 'adapters.txt' file contains common universal adapters for Illumina platforms, but can be modified to fit your specific protocol. Once the reads are trimmed, any reads that are shorter than 25 bp are removed. Single-end sequences are then output as '_trimmed.fastq(.gz)' files. If using Paired-end samples, Trimmomatic checks pairing of forward and reverse reads and outputs '_paired' and '_unpaired' files for each input '.fastq(.gz)' file. See [project website](http://www.usadellab.org/cms/?page=trimmomatic) for more details.

- **Alignment**:

    The '_trimmed' and/or '_paired' files are then aligned to the parent genome using Bowtie2. Bowtie2 indexes are generated from the '.fasta' file provided by the 'genome [-g]' argument if not already available. To ensure that all possible alignments are identified, the --very-sensitive option is used to specify a high-sensitivity mode. A single SAM file is output for each aligned sample. See [project website](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for more details.

- **Deduplication**:

    PCR duplicates are removed from the SAM formatted '_aligned.sam' genome alignments. Deduplication metrics are output in a text file. See [project website](https://broadinstitute.github.io/picard/) for more details.

- **Filtering**:

    Low-quality reads with Phred score less than 30 are removed from the deduplicated '_dedup.sam' files using Samtools view. The output BAM format '_dedup.bam' file only includes properly paired aligned reads using the flags '-f 0x02' and '-F 0x04'. See [project website](https://www.htslib.org/doc/samtools.html) for more details.

- **Sorting**:

    Once deduplicated and filtered, the '_dedup.bam' files are sorted by coordinate using Samtools sort to make mapping faster. See [project website](https://www.htslib.org/doc/samtools.html) for more details.

- **Mapping**:

    To get the depth of coverage at each nucleotide, the '_sorted.bam' files are mapped to the genome using Samtools depth. The output BED file is a tab-delimited file with three columns: 'chromosome', 'position' and 'reads'. See [project website](https://www.htslib.org/doc/samtools.html) for more details.

- **Peak calling**:

    Peak calling is performed using macs2 callpeak. 
