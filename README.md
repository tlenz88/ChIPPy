# ChIPPy
This is a complete pipeline for ChIP-seq data analysis. Pipeline can be started from any step, but specific files are required for the desired step.

It is important to note that ChIPPy and any supplementary tools are merely templates to provide fast and simple ChIP-seq analysis and may not work for all samples or experimental designs. If new to ChIP-seq analysis, please learn about each of the individual tools used in this package and cite their respective publications.

## How do I organize my data?

Put all input files in a single directory with one folder per sample. The sequences can be either single or paired-end and can be gzipped, though this is not a requirement. Paired-end files need to have '_R1' and '_R2' within the file names to specify forward and reverse reads.

```
    + PATH_TO_INPUT
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

## What tools do I need?

Although ChIPPy automatically downloads most of the necessary software/tools, depending on the system environment the user may need to install some things manually. The only software that requires manual installation is [**Trimmomatic**](http://www.usadellab.org/cms/?page=trimmomatic) due to compatibility issues with the java version specified in the Trimmomatic build.xml file. Here is a list of all the necessary tools to run the complete pipeline:

- [Python3](https://www.python.org/downloads/)
- [pip](https://pip.pypa.io/en/stable/installation/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Picard](https://broadinstitute.github.io/picard/)
- [Samtools](https://www.htslib.org/doc/samtools.html)

### Conda environment installation

To make installation of ChIPPy dependencies easier, users can create a conda environment using the provided YAML file. Any dependencies for other utilities within this package will also be installed.

1. Download and install [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [anaconda](https://www.anaconda.com/download).

2. Create conda environment: ```conda env create -f environment.yml```

3. Activate environment: ```conda activate ChIPPy```

## What else do I need?

Any genome can be used if a properly formatted FASTA file is provided. Any additional files, such as bowtie2 indexes, will be created automatically if not already present.

The user may want to download [IGV](https://igv.org/) because a WIG file will be generated for easier visualization in the IGV genome browser.

Additional scripts used to automate plotting are available to the user in the ```utils``` folder if a different pipeline or tools were used to perform ChIP-seq analysis.

## How do I run ChIPPy?

Brief description of input arguments via help message:

```
    ChIPPy.sh --help
    Usage: ChIPPy.sh -i INPUT -g GENOME [-o OUTPUT] [-s STEP] [-q QUALITY] [-t THREADS] [-r]

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
      -t|--threads THREADS     : Processor threads.
      -r|--remove REMOVE       : Remove intermediate files.
      -h|--help HELP           : Show help message.
    ----------------------------------------------------------------
```

The following is a more detailed description of input arguments:

- **Required arguments**:

    --input [-i]: Directory containing a folder for each sample. Samples can be single or paired-end, but paired-end sequences need to contain '_R1' and '_R2' within the file names to determine forward and reverse reads. If starting from a step in the pipeline after alignment the input files are not required to contain '_R1' or '_R2'.

    --genome [-g]: Directory containing genome files for the organism of interest. A FASTA file should be present in the directory so that any additional required genome files can be automatically created if missing. All genome files should have the same basename.

- **Optional arguments**:

    --output [-o]: Output directory. If no output is given, the output files will be saved to the input folders.

    --step [-s]: Desired step at which to start the pipeline. If no step is entered, the entire pipeline will be run.

    --quality [-q]: Integer indicating Phred quality score for trimming low-quality bases at the ends of reads during read pairing and for removing low-quality reads during filtering. Phred score will also be used to approximate filtering settings during alignment.

    --threads [-t]: Integer indicating number of processor threads to use for tools that allow multithreading. If no value is given, the number of available threads will be determined automatically and will use half of available threads on the users system.

    --remove [-r]: Removes intermediate files when no longer needed by pipeline to save hard drive space.

## What are the individual steps within the pipeline?

- **Quality check**:

    FastQC is used to check sequence quality prior to running analysis pipeline. Statistics such as per base sequence quality, per base sequence content, sequence duplication levels and adapter content are calculated from input FASTQ files. See the [project website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for more details.

- **Trimming**:

    Trimmomatic is used to perform initial trimming. Adapter sequences in the provided 'adapters.txt' file and bases below the quality threshold set by the ```quality [-q]``` argument (default = 30) are removed from the ends of reads. The 'adapters.txt' file contains common universal adapters for Illumina platforms, but can be modified to fit a specific protocol. Once the reads are trimmed, any reads that are shorter than 25 bp are removed. Single-end sequences are then output as trimmed FASTQ files. If using paired-end samples, Trimmomatic checks pairing of forward and reverse reads and outputs paired and unpaired files for each input FASTQ file. See the [project website](http://www.usadellab.org/cms/?page=trimmomatic) for more details.

- **Alignment**:

    The trimmed and/or paired files are then aligned to the parent genome using Bowtie2. Bowtie2 indexes are generated from the FASTA file provided by the ```genome [-g]``` argument if not already available. To ensure that all possible alignments are identified, the --very-sensitive argument is used to specify the high-sensitivity mode. A single SAM file is output for each aligned sample. See the [project website](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for more details.

- **Deduplication**:

    PCR duplicates are removed from the SAM files using ```picard markDuplicates```. Deduplication metrics are output in a text file. See the [project website](https://broadinstitute.github.io/picard/) for more details.

- **Filtering**:

    Low-quality reads with Phred score less than 30 are removed from the deduplicated SAM files using ```samtools view```. The filtered BAM files only includes properly paired aligned reads using the flags ```-f 0x02 -F 0x04```. See the [project website](https://www.htslib.org/doc/samtools.html) for more details.

- **Sorting**:

    Once deduplicated and filtered, the filtered BAM files are sorted by coordinate using ```samtools sort``` to make mapping faster. See the [project website](https://www.htslib.org/doc/samtools.html) for more details.

- **Mapping**:

    To get the depth of coverage at each nucleotide, the sorted BAM files are mapped to the genome using ```samtools depth```. The output BED file is a tab-delimited file with three columns: chromosome name, position and reads. See the [project website](https://www.htslib.org/doc/samtools.html) for more details.


# ChIPPeaks

ChIPPeaks is a supplementary tool that performs peak calling and differential peak calling. It is meant to be used after running ChIPPy on a set of samples, but can be used on any dataset if the proper inputs are provided.

Peak calling is performed using ```macs2 callpeak```. The ```q-value [-q]``` is set at 0.05 and ```genome size [-g]``` is determined automatically. If the data is paired-end, the format of the input will be set to ```-f BEDPE``` so that the insert size of pairs is used to build fragment pileup. The choice of peak calling algorithm--broad or narrow--is also determined automatically based on the 'Factor' column of the metadata file. Broad peak calling is performed for histone modifications, whereas narrow peak calling is for transcription factors. See the [MACS github repo](https://github.com/macs3-project/MACS) for more details.

The DiffBind R package is used for differential peak calling analysis. The metadata file is used to dictate experimental design and will determine what comparisons are made. See the [DiffBind vignette](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) for more information.

## How should I set up the metadata?

The metadata file should contain a minimum of 7 columns used to describe the format of the user's ChIP-seq experiment, including sample names, relationships between samples, replicate numbers, and file paths. Although it isn't required for the metadata to be formatted exactly as shown in the ```example_metadata.txt``` file, it's recommended that any additional columns be added at the end (to the right) of the existing columns. The 'Condition' column is used for differential peak calling.

The following are brief descriptions of the template columns:

    --SampleID: Short string unique to each target sample.

    --Factor: Antibody used for immunoprecipitation of each target sample.

    --Condition: String used to differentiate sample conditions, i.e. knockout (KO) vs control (CTRL)

    --Replicate: Integer indicating replicate number.

    --bamReads: Path to coordinate sorted BAM files for target samples.

    --ControlID: Short strings unique to each control sample.

    --bamControl: Path to coordinate sorted BAM files for control samples.

## ChIPPeaks.sh

Brief description of input arguments via help message:

```
    ChIPPeaks.sh --help
    Usage: ChIPPeaks.sh -m METADATA -c CONTROL -g GENOME [-h]

    --------------------------------------------------------------
     Input arguments:
      -m|--metadata METADATA : ChIP-seq sample metadata file.
      -c|--control CONTROL   : Control 'Condition'.
      -g|--genome GENOME     : Directory containing genome files.
      -h|--help HELP         : Show help message.
    --------------------------------------------------------------
```

The following is a more detailed description of input arguments:

- **Input arguments**:

    --metadata [-m]: Path to the tab-delimited text file containing ChIP-seq experiment metadata. Columns of metadata file should have at least seven columns: SampleID, Factor, Condition, Replicate, bamReads, ControlID and bamControl. The metadata file will be automatically modified to include two additional columns--Peaks and PeakCaller--after peak calling. Use the ```example_metadata.txt``` file in the ```examples``` directory as a template. See the [DiffBind vignette](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf) and the following section of this README for more information on setting up the metadata file.

    --control[-c]: Sting indicating the control samples in the experiment. This string should be under the 'Condition' column in the metadata file.

    --genome [-g]: Directory containing genome files for the organism of interest. A FASTA file should be present in the directory so that any additional required genome files can be automatically created if missing. All genome files should have the same basename.
