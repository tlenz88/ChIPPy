# ChIPPy Utilities
ChIPPy uses several python scripts to generate input files or automate plotting of ChIP-seq data. The plotting tools can easily be used independently if a different pipeline or tools were used to perform ChIP-seq analysis.

## plot_chromosome_coverage.py

Plots genomewide coverage for ChIP-seq data. Track lengths are normalized with respect to the longest chromosome--i.e. the longest chromosome will fill the width of the figure and all other chromosomes are proportionally plotted against it.

**Usage**: ```plot_chromosome_coverage.py  -b sample1.bed sample2.bed sample3.bed```

- **Required input arguments**:

    --bed [-b]: BED file(s) containing genome-wide per-base coverage. To generage this file from a BAM file: ```samtools depth -a -o input.bed input.bam```

- **Optional input arguments**:

    --gff [-g]: Standard GFF file for genome/organism of interest.

    --samples [-s]: List of sample names to annotate the y-axes of barplots. List should be the same length as the list of BED files and plots will be annotated in the same order as BED files. If no sample names are given, the BED file names are used to annotate plots.

    --output [-o]: Directory and name for the output PDF file. If no output is given, the output file is 'chromosome_coverage.pdf' and saved in the directory of the first BED input.

    --resolution [-r]: Integer indicating the binning resolution (default = 1000).

    --centromeres [-c]: Tab-delimited file containing centromere coordinates for each chromosome. File should have three columns: chromosome names, start coordinates and end coordinates. The values in the 'chromosome names' column should match the chromosome names in the BED files.

    --gene_list [-l]: List of genes to plot, separated by spaces, commas, newlines, or tabs. Genes will be extracted from the provided GFF file. If no list is provided, all genes will be annotated on plots.

    --normalization [-n]: If argument is provided, data is counts-per-million (CPM) normalized prior to plotting. This allows for direct comparison of samples regardless of sequencing depth.

    --ymax [-y]--If argument is provided, all chromosome plots will use the maximum y-value for the entire binned dataset. If not provided, each chromosome will use the maximum value for that chromosome.

## plot_gene_coverage.py

**Usage**: ```plot_gene_coverage.py -b sample1.bed sample2.bed sample3.bed -g genome.gff```

- **Required input arguments**:

    --bed [-b]: BED file(s) containing genome-wide per-base coverage. To generage this file from a BAM file: ```samtools depth -a -o input.bed input.bam```

    --gff [-g]: Standard GFF file for genome/organism of interest.

- **Optional input arguments**:

    --samples [-s]: List of sample names to annotate the y-axes of barplots. List should be the same length as the list of BED files and plots will be annotated in the same order as BED files. If no sample names are given, the BED file names are used to annotate plots.

    --output [-o]: Directory and name for the output PDF file. If no output is given, the output file is 'chromosome_coverage.pdf' and saved in the directory of the first BED input.

    --resolution [-r]: Integer indicating the binning resolution (default = 10).

    --gene_list [-l]: List of genes to plot, separated by spaces, commas, newlines, or tabs. Genes will be extracted from the provided GFF file. If no list is provided, all genes will be annotated on plots.

    --normalization [-n]: If argument is provided, data is counts-per-million (CPM) normalized prior to plotting. This allows for direct comparison of samples regardless of sequencing depth.

    --distance [-d]: Integer indicating length of 5' and 3' regions to plot. By default, only the coding region for the plotted genes is used.

    --ymax [-y]--If argument is provided, all chromosome plots will use the maximum y-value for the entire binned dataset. If not provided, each chromosome will use the maximum value for that chromosome.

## ChIPPeaks.sh

Peak calling is performed using ```macs2 callpeak```. The ```q-value [-q]``` is set at 0.05 and ```genome size [-g]``` is determined automatically. If the data is paired-end, the format of the input will be set to ```-f BEDPE``` so that the insert size of pairs is used to build fragment pileup. The choice of peak calling algorithm--broad or narrow--is also determined automatically based on the input bam file being analyzed. Broad peak calling is performed for histone modifications, whereas narrow peak calling is for transcription factors. If no peaks are called using the default ```--mfold``` parameter, ```macs2 callpeak``` will be repeatedly run using a decreasing lower limit for model building until enough peaks are found to build the shifting model or until the lower limit reaches 1. See the [github repo](https://github.com/macs3-project/MACS) for more details.