
### This is a complete pipeline for ChIP-seq data analysis. Required packages are installed automatically. Pipeline can be started from any step, but specific files are required for the desired step.

## Organize your data

Put all input files in a single directory with one folder per sample, such as:

```
    + PATH_TO_DATA
        + sample1
            ++ sample1_R1.fastq.gz
            ++ sample1_R2.fastq.gz
        + sample2
            ++ sample2_R1.fastq.gz
            ++ sample2_R2.fastq.gz
```

## How do I use it?

Brief description of input arguments via help message:

```
    ChIPpipe.sh --help
    usage : ChIPpipe.sh -i INPUT -g GENOME [-o OUTPUT] [-s STEP] [-q QUALITY] [-t THREADS] [-r]

    -------------------------------------------------------------
    Required inputs:
     -i|--input  INPUT    : Input data folder.
     -g|--genome GENOME   : Path to genome files.

    Optional inputs:
     -o|--output OUTPUT   : Output folder.
     -s|--step  STEP      : Choose starting step.
            quality_check : Initial quality check.
                 trimming : Adapter trimming.
                alignment : Read alignment.
            deduplication : Remove PCR duplicates.
                filtering : Filtering low-quality reads.
                  sorting : Sorting reads by coordinate.
                  mapping : Mapping reads to each base pair
     -q|--quality QUALITY : Phred quality score for filtering.
     -t|--threads THREADS : Processor threads.
     -r|--remove REMOVE   : Remove intermediate files.
    -------------------------------------------------------------
```

The following is a more detailed description of input arguments:

- Required arguments:
    <p>--input [-i]:  Directory containing a folder for each sample. Samples can be single or paired-end, but paired-end sequences need to contain '_R1' and '_R2' within the file names to determine forward and reverse reads. If starting from a step in the pipeline after alignment the input files are not required to contain '_R1' or '_R2'.</p>
    <p>--genome [-g]: Directory containing genome files for the organism of interest. A '.fasta' file should be present in the directory so that any additional required genome files can be automatically created if missing. All genome files should have the same basename.</p>

- Optional arguments:
    <p>--output [-o]:  Output directory. If no output is given, the output files will be saved to the input folders.</p>
    <p>--step [-s]:    Desired step at which to start the pipeline. If no step is entered, the entire pipeline will be run.</p>
    <p>--quality [-q]: Integer indicating Phred quality score for trimming low-quality bases at the ends of reads during read pairing and for removing low-quality reads during filtering. Phred score will also be used to approximate filtering settings during alignment.</p>
    <p>--threads [-t]: Integer indicating number of processor threads to use for tools that allow multithreading. If no value is given, the number of available threads will be determined automatically and will use half of available threads on your system.</p>
    <p>--remove [-r]:  Removes intermediate files when no longer needed by pipeline to save hard drive space.
