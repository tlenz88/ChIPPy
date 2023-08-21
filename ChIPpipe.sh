#!/bin/bash

## Complete pipeline for ChIP-seq data analysis.
## Created June 13, 2022
## Updated August 17, 2023
## Author(s): Todd Lenz
## Contact: tlenz001@ucr.edu

function help {
    echo "ChIPpipe.sh --help"
    echo "usage : ChIPpipe.sh -i INPUT -g GENOME [-o OUTPUT] [-s STEP] [-q QUALITY] [-a1 ANTIBODY1] [-a2 ANTIBODY2] [-t THREADS] [-r]"
    echo
    echo "----------------------------------------------------------------"
    echo "Required inputs:"
    echo "  -i|--input  INPUT        : Input data folder."
    echo "  -g|--genome GENOME       : Path to genome files."
    echo
    echo "Optional inputs:"
    echo "  -o|--output OUTPUT       : Output folder."
    echo "  -s|--step STEP           : Choose starting step."
    echo "         quality_check     : Initial quality check."
    echo "              trimming     : Adapter trimming."
    echo "             alignment     : Read alignment."
    echo "         deduplication     : Remove PCR duplicates."
    echo "             filtering     : Filtering low-quality reads."
    echo "               sorting     : Sorting reads by coordinate."
    echo "               mapping     : Mapping reads to each base pair"
    echo "  -q|--quality QUALITY     : Phred quality score for filtering."
    echo " -a1|--antibody1 ANTIBODY1 : Target antibody."
    echo " -a2|--antibody2 ANTIBODY2 : Control antibody (IGG) or input."
    echo "  -t|--threads THREADS     : Processor threads."
    echo "  -r|--remove REMOVE       : Remove intermediate files."
    echo "  -h|--help HELP           : Show help message."
    echo "----------------------------------------------------------------"
    exit 0;
}


#############################
## Define input arguments. ##
#############################
for arg in "$@"; do
    case "$arg" in
        "--input") set -- "$@" "-i" ;;
        "--genome") set -- "$@" "-g" ;;
        "--output") set -- "$@" "-o" ;;
        "--step") set -- "$@" "-s" ;;
        "--quality") set -- "$@" "-q" ;;
        "--antibody1") set -- "$@" "-a1" ;;
        "--antibody2") set -- "$@" "-a2" ;;
        "--threads") set -- "$@" "-t" ;;
        "--remove") set -- "$@" "-r" ;;
        "--help") set -- "$@" "-h" ;;
        *) set -- "$@" "$arg"
    esac
done

pipedir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
SCRIPTS="$pipedir"/scripts

while getopts ":i:g:o:s:q:a1:a2:t:r:h:" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        g) GENOME="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        s) STEP="$OPTARG";;
        q) QUALITY="$OPTARG";;
        a1) ANTIBODY1="$OPTARG";;
        a2) ANTIBODY2="$OPTARG";;
        t) THREADS="$OPTARG";;
        r) REMOVE=true;;
        h) help;;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done


#############################
## Define input and output ##
#############################
if [[ -e "$OUTPUT" && "$STEP" == "" ]]; then
    echo "$OUTPUT folder alreads exists. Do you want to overwrite it? (y/n)"
    read -r ans
    if [[ "$ans" == "y" ]]; then
        rm -rf "$OUTPUT"
        mkdir "$OUTPUT"
        for i in "$INPUT"/*; do
            mkdir "$OUTPUT"/"$(basename "$i")"
        done
    fi
elif [[ "$OUTPUT" != "" && "$STEP" == "" ]]; then
    mkdir "$OUTPUT"
    for i in "$INPUT"/*; do
        mkdir "$OUTPUT"/"$(basename "$i")"
    done
elif [[ "$OUTPUT" == "" && "$STEP" == "" ]]; then
    echo "No output folder selected. Files will be saved to input."
    OUTPUT=$INPUT
fi


######################################################
## Check Python installation and required packages. ##
######################################################
for p in "python3" "python3-pip" "fastqc" "bowtie2" "samtools"; do
    if ! command -v "$p" &> /dev/null; then
        echo "Installing newest version of $p."
        if command -v apt &> /dev/null; then
            sudo apt update
            sudo apt upgrade
            sudo apt install -y "$p"
        elif command -v yum &> /dev/null; then
            sudo yum install -y "$p"
        elif command -v dnf &> /dev/null; then
            sudo dnf update
            sudo dnf install -y "$p"
        elif command -v zypper &> /dev/null; then
            sudo zypper refresh
            sudo zypper update
            sudo zypper --non-interactive install "$p"
        elif command -v pacman &> /dev/null; then
            sudo pacman -Syu
            sudo pacman install --noconfirm "$p"
        elif command -v brew &> /dev/null; then
            brew update
            brew upgrade
            brew install -y "$p"
        else
            echo "Can't determine package manager. Install $p manually."
            exit 1
        fi
    fi
done

if ! command -v macs2 &> /dev/null; then
    echo "Installing macs2."
    pip3 install macs2
fi

for pkg in "biopython" "matplotlib" "numpy" "pandas"; do
    if ! python3 -c "import $pkg" &> /dev/null; then
        echo "Installing $pkg."
        pip3 install "$pkg"
    fi
done

if [[ $(find -L "$SCRIPTS" -mindepth 1 -maxdepth 4 -name "picard.jar" | wc -l) == 0 ]]; then
    if [ ! -d "$SCRIPTS"/picard ]; then
        echo "Installing picard."
        git clone https://github.com/broadinstitute/picard.git "$SCRIPTS"/picard
    fi
    "$SCRIPTS"/picard/gradlew -p "$SCRIPTS"/picard shadowJar
fi


###########################################
## Check for all necessary genome files. ##
###########################################
if [[ $GENOME != "" ]]; then
    GENOME="$(readlink -f "$GENOME")"
    fa=$(find -L "$GENOME" -mindepth 1 -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \))
    bname="${fa%%.*}"
    if [[ ! -e "$bname.1.bt2" ]]; then
        echo "Creating Bowtie2 indexes."
        bowtie2-build --threads "$THREADS" "$fa" "$bname"
    else
        echo "Bowtie2 indexes detected."
    fi
    if [[ ! -e "$bname.chrom.sizes" ]]; then
        echo "Creating chromosome sizes file."
        "$SCRIPTS"/create_sizes.py "$fa"
    else
        echo "Chromosome sizes file detected."
    fi
    if [[ ! -e "$fa.fai" ]]; then
        echo "Creating FASTA index."
        samtools faidx "$fa"
    else
        echo "FASTA index detected."
    fi
    if [[ ! -e "$bname.gff" ]]; then
        echo "GFF file not detected. Some plots will not be generated."
    else
        echo "GFF file detected."
    fi
else
    echo "Error: No input genome (-g) detected."
    exit 1
fi


#################################################
## Define number of processing threads to use. ##
#################################################
if [[ -z $THREADS ]]; then
    THREADS=$(($(nproc) / 2))
fi


####################################################
## Check starting step and required files option. ##
####################################################
STEP="${STEP/^ //}"
ALL_STEPS=("quality_check" "trimming" "alignment" "deduplication" "filtering" "sorting" "mapping" "peak_calling")
QUALITY_CHECK="quality_check"
TRIMMING="trimming"
ALIGNMENT="alignment"
DEDUPLICATION="deduplication"
FILTERING="filtering"
SORTING="sorting"
MAPPING="mapping"
PEAK_CALLING="peak_calling"

if [[ $STEP != "" ]]; then
    case $STEP in
        "quality_check") FILTERED_STEPS=("${ALL_STEPS[@]}");;
        "trimming") FILTERED_STEPS=("${ALL_STEPS[@]:1}");;
        "alignment") FILTERED_STEPS=("${ALL_STEPS[@]:2}");;
        "deduplication") FILTERED_STEPS=("${ALL_STEPS[@]:3}");;
        "filtering") FILTERED_STEPS=("${ALL_STEPS[@]:4}");;
        "sorting") FILTERED_STEPS=("${ALL_STEPS[@]:5}");;
        "mapping") FILTERED_STEPS=("${ALL_STEPS[@]:6}");;
        "peak_calling") FILTERED_STEPS=("${ALL_STEPS[@]:7}");;
        *) echo "Invalid step: $STEP"; exit 1;;
    esac
else
    FILTERED_STEPS=("${ALL_STEPS[@]}")
fi

if [[ "$STEP" != "" ]]; then
    for s in $STEP; do
        echo
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
        for i in "${QUALITY_CHECK[@]}"; do
            if [[ "$i" == "$s" ]]; then
                fq=$(find -L "$INPUT" -mindepth 2 -maxdepth 2  -name "*.fastq*" -o -name "*.fastq*" | wc -l)
                if [[ "$fq" == 0 ]]; then
                    echo "Error: No '.fastq(.gz)' files detected."
                    exit 1
                elif [[ "$fq" == 2 ]]; then
                    fq_r1=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R1*" | wc -l)
                    fq_r2=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R2*" | wc -l)
                    if [[ "$fq_r1" == 0 || "$fq_r2" == 0 || "$fq_r1" != "$fq_r2" ]]; then
                        echo "Error: Paired-end input files must contain '_R1' and '_R2' in the file names."
                        exit 1
                    fi
                fi
            fi
        done
        for i in "${TRIMMING[@]}"; do
            if [[ "$i" == "$s" ]]; then
                fq=$(find -L "$INPUT" -mindepth 2 -maxdepth 2  -name "*.fastq*" -o -name "*.fastq*" | wc -l)
                if [[ "$fq" == 0 ]]; then
                    echo "Error: No '.fastq(.gz)' files detected."
                    exit 1
                elif [[ "$fq" == 2 ]]; then
                    fq_r1=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R1*" | wc -l)
                    fq_r2=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R2*" | wc -l)
                    if [[ "$fq_r1" == 0 || "$fq_r2" == 0 || "$fq_r1" != "$fq_r2" ]]; then
                        echo "Error: Paired-end input files must contain '_R1' and '_R2' in the file names."
                        exit 1
                    fi
                fi
            fi
        done
        for i in "${ALIGNMENT[@]}"; do
            if [[ "$i" == "$s" ]]; then
                tfq=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 -name "*_paired.fastq*" -o -name "*_paired.fq*" -o -name "*_trimmed.fastq*" -o -name "*_trimmed.fq*" | wc -l)
                if [[ "$tfq" == 0 ]]; then
                    echo "Error: Trimmed fastq files are required to start the pipeline at the alignment step."
                    exit 1
                elif [[ "$tfq" == 1 ]]; then
                    tfq_r1=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*_trimmed.fastq*" -o -name "*_trimmed.fq*" \) -and -name "*_R1*" | wc -l)
                    if [[ "$tfq_r1" == 0 ]]; then
                        echo "Error: Trimmed fastq files are required to start the pipeline at the alignment step."
                        echo "Ensure that single-end input files contain '_trimmed.fastq' and '_R1' in the file name."
                        exit 1
                    fi
                elif [[ "$tfq" == 2 ]]; then
                    tfq_r1=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*_paired.fastq*" -o -name "*_paired.fq*" \) -and -name "*_R1*" | wc -l)
                    tfq_r2=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*_paired.fastq*" -o -name "*_paired.fq*" \) -and -name "*_R2*" | wc -l)
                    if [[ "$tfq_r1" == 0 || "$tfq_r2" == 0 || "$tfq_r1" != "$tfq_r2" ]]; then
                        echo "Error: Trimmed fastq files are required to start the pipeline at the alignment step."
                        echo "Ensure that paired-end input files contain '_paired.fastq' and '_R1'/'_R2' in the file name."
                        exit 1
                    fi
                fi
            fi
        done
        for i in "${DEDUPLICATION[@]}"; do
            if [[ "$i" == "$s" ]]; then
                asam=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 -name "*_aligned.sam" | wc -l)
                if [[ "$asam" == 0 ]]; then
                    echo "Error: Aligned SAM files are required to start the pipeline at the deduplication step."
                    exit 1
                fi
            fi
        done
        for i in "${FILTERING[@]}"; do
            if [[ "$i" == "$s" ]]; then 
                dsam=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 -name "*_dedup.sam" | wc -l)
                if [[ "$dsam" == 0 ]]; then
                    echo "Error: Deduplicated SAM files are required to start the pipeline at the filtering step."
                    exit 1
                fi
            fi
        done
        for i in "${SORTING[@]}"; do
            if [[ "$i" == "$s" ]]; then 
                fbam=$(find -L "$INPUT" -mindepth 2 -maxdepth 4 -name "*_dedup.bam" | wc -l)
                if [[ "$fbam" == 0 ]]; then
                    echo "Error: Filtered BAM files are required to start the pipeline at the sorting step."
                    exit 1
                fi
            fi
        done
        for i in "${MAPPING[@]}"; do
            if [[ "$i" == "$s" ]]; then 
                sbam=$(find -L "$INPUT" -mindepth 2 -maxdepth 4 -name "*_sorted.bam" | wc -l)
                if [[ "$sbam" == 0 ]]; then
                    echo "Error: Sorted BAM files are required to start the pipeline at the mapping step."
                    exit 1
                fi
            fi
        done
        for i in "${PEAK_CALLING[@]}"; do
            if [[ "$i" == "$s" ]]; then 
                sbam=$(find -L "$INPUT" -mindepth 2 -maxdepth 4 -name "*_sorted.bam" | wc -l)
                if [[ "$sbam" == 0 ]]; then
                    echo "Error: Sorted BAM files are required to start the pipeline at the peak_calling step."
                    exit 1
                fi
            fi
        done
    done
else
    echo
    echo "Starting ChIP-seq analysis pipeline."
    echo
fi


###################################################
## List of functions/steps in ChIP-seq pipeline. ##
###################################################
function quality_check() {
    echo "Performing FastQC quality check."
    for i in "$INPUT"/*; do
        fq=$(find -L "$i" -mindepth 1 -maxdepth 1  -name "*.fastq*" -o -name "*.fastq*" | wc -l)
        if [[ "$fq" == 1 ]]; then
            fq_r1=$(find -L "$i" -mindepth 1 -maxdepth 1  -name "*.fastq*" -o -name "*.fastq*" )
            fastqc -t "$THREADS" "$fq_r1"
            if [[ "$OUTPUT" != "$INPUT" ]]; then
                mv "${fq_r1%%.*}"_fastqc.html "$OUTPUT/$(basename "$(dirname "$fq_r1")")/$(basename "${fq_r1%%.*}")_fastqc.html"
            fi
            rm "${fq_r1%%.*}"_fastqc.zip
        fi
        if [[ "$fq" == 2 ]]; then
            fq_r1=$(find -L "$i" -mindepth 1 -maxdepth 1 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R1*")
            fq_r2=$(find -L "$i" -mindepth 1 -maxdepth 1 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R2*")
            fastqc -t "$THREADS" "$fq_r1"
            fastqc -t "$THREADS" "$fq_r2"
            if [[ "$OUTPUT" != "$INPUT" ]]; then
                mv "${fq_r1%%.*}"_fastqc.html "$OUTPUT/$(basename "$(dirname "$fq_r1")")/$(basename "${fq_r1%%.*}")_fastqc.html"
                mv "${fq_r2%%.*}"_fastqc.html "$OUTPUT/$(basename "$(dirname "$fq_r2")")/$(basename "${fq_r2%%.*}")_fastqc.html"
            fi
            rm "${fq_r1%%.*}"_fastqc.zip
            rm "${fq_r2%%.*}"_fastqc.zip
        elif [[ "$fq" == 0 ]]; then
            echo "Error: Input files must be '.fastq(.gz)'."
            exit 1
        fi
    done
    echo
}

function trimming() {
    if [[ "$STEP" == "trimming" ]]; then
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
    echo "Trimming and pairing R1 and R2 reads."
    for i in "$INPUT"/*; do
        basename "$i"
        fq=$(find -L "$i" -mindepth 1 -maxdepth 1  -name "*.fastq*" -o -name "*.fastq*" | wc -l)
        tmjar=$(find -L "$SCRIPTS" -mindepth 1 -maxdepth 3 -name "*trimmomatic*")
        if [[ "$fq" == 1 ]]; then
            fq_r1=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*.fastq*" -o -name "*.fq*" )
            out_r1="$OUTPUT"/"$(basename "$(dirname "$fq_r1")")"/"$(basename "${fq_r1%%.*}")"
            java -jar "$tmjar" SE -threads "$THREADS" "$fq_r1" "${out_r1}_trimmed.fastq.gz" LEADING:"$QUALITY" TRAILING:"$QUALITY" MINLEN:25 ILLUMINACLIP:"$GENOME"/Illumina_CD_index_adapters.txt:2:30:10
        elif [[ "$fq" == 2 ]]; then
            fq_r1=$(find -L "$i" -mindepth 1 -maxdepth 1 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R1*")
            fq_r2=$(find -L "$i" -mindepth 1 -maxdepth 1 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R2*")
            out_r1="$OUTPUT"/"$(basename "$(dirname "$fq_r1")")"/"$(basename "${fq_r1%%.*}")"
            out_r2="$OUTPUT"/"$(basename "$(dirname "$fq_r2")")"/"$(basename "${fq_r2%%.*}")"
            java -jar "$tmjar" PE -threads "$THREADS" "$fq_r1" "$fq_r2" "${out_r1}_paired.fastq.gz" "${out_r1}_unpaired.fastq.gz" "${out_r2}_paired.fastq.gz" "${out_r2}_unpaired.fastq.gz" LEADING:"$QUALITY" TRAILING:"$QUALITY" MINLEN:25 ILLUMINACLIP:"$GENOME"/Illumina_CD_index_adapters.txt:2:30:10
        elif [[ "$fq" == 0 ]]; then
            echo "Error: Input files must be '.fastq(.gz)'."
            exit 1
        fi
    done
    echo
}

function alignment() {
    if [[ "$STEP" == "alignment" ]]; then
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
        DIR="$INPUT"
    else
        DIR="$OUTPUT"
    fi
    echo "Aligning reads to the genome."
    for i in "$DIR"/*; do
        basename "$i"
        tfq=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_trimmed.fastq*" -o -name "*_trimmed.fq*" -o -name "*_paired.fastq*" -o -name "*_paired.fq*"  | wc -l)
        fa="$(find -L "$(readlink -f "$GENOME")" -mindepth 1 -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \))"
        bidx="${fa%%.*}"
        if [[ $tfq == 1 ]]; then
            tfq_r1=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_trimmed.fastq*" -o -name "*_trimmed.fq*" )
            out="$OUTPUT/$(basename "$(dirname "$(readlink -f "$tfq_r1")")")/$(basename "${tfq_r1%_*}")"
            bowtie2 -x "$bidx" -1 "$tfq_r1" -S "${out}_aligned.sam" -p "$THREADS" --very-sensitive
            if [[ "$REMOVE" == true ]]; then
                rm "$i"/*_trimmed.fastq*
            fi
        elif [[ $tfq == 2 ]]; then
            tfq_r1=$(find -L "$i" -mindepth 1 -maxdepth 1 \( -name "*_paired.fastq*" -o -name "*_paired.fq*" \) -and -name "*_R1*")
            tfq_r2=$(find -L "$i" -mindepth 1 -maxdepth 1 \( -name "*_paired.fastq*" -o -name "*_paired.fq*" \) -and -name "*_R2*")
            out="$OUTPUT/$(basename "$(dirname "$(readlink -f "$tfq_r1")")")/$(basename "${tfq_r1%_*}")"
            asam="${out/_R1/}"
            bowtie2 -x "$bidx" -1 "$tfq_r1" -2 "$tfq_r2" -S "${asam}_aligned.sam" -p "$THREADS" --very-sensitive
            if [[ "$REMOVE" == true ]]; then
                rm "$i"/*_paired.fastq* "$i"/*_unpaired.fastq*
            fi
        elif [[ "$tfq" == 0 ]]; then
            echo "Error: Trimmed '.fastq(.gz)' files are required for the alignment step."
            exit 1
        fi
    done
    echo
}

function deduplication() {
    if [[ "$STEP" == "deduplication" ]]; then
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
        DIR="$INPUT"
    else
        DIR="$OUTPUT"
    fi
    echo "Removing PCR duplicates."
    for i in "$DIR"/*; do
        basename "$i"
        sam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_aligned.sam" | wc -l)
        ptjar=$(find -L "$SCRIPTS" -mindepth 1 -maxdepth 4 -name "picard.jar")
        if [[ "$sam" == 1 ]]; then
            asam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_aligned.sam")
            dsam="$OUTPUT/$(basename "$(dirname "$(readlink -f "$asam")")")/$(basename "${asam%_*}")"
            java -jar "$ptjar" MarkDuplicates -I "$asam" -O "${dsam}_dedup.sam" -M "${dsam}_metrics.txt" -ASO queryname
            if [[ "$REMOVE" == true ]]; then
                rm "$asam"
            fi
        elif [[ "$sam" == 0 ]]; then
            echo "Error: SAM files are required for the deduplication step."
            exit 1
        fi
    done
    echo
}

function filtering() {
    if [[ "$STEP" == "filtering" ]]; then
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
        DIR="$INPUT"
    else
        DIR="$OUTPUT"
    fi
    echo "Filtering low quality reads."
    for i in "$DIR"/*; do
        basename "$i"
        sam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_dedup.sam" | wc -l)
        if [[ "$sam" == 1 ]]; then
            dsam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_dedup.sam")
            dbam="$OUTPUT/$(basename "$(dirname "$(readlink -f "$dsam")")")/$(basename "${dsam%%.*}")"
            samtools view -@ "$THREADS" -q "$QUALITY" -f 0X02 -F 0X04 -b "$dsam" -o "${dbam}.bam"
            if [[ "$REMOVE" == true ]]; then
                rm "$dsam"
            fi
        elif [[ "$sam" == 0 ]]; then
            echo "Error: Deduplicated SAM files are required for the filtering step."
            exit 1
        fi
    done
    echo
}

function sorting() {
    if [[ "$STEP" == "sorting" ]]; then
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
        DIR="$INPUT"
    else
        DIR="$OUTPUT"
    fi
    echo "Sorting reads"
    for i in "$DIR"/*; do
        basename "$i"
        bam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_dedup.bam" | wc -l)
        if [[ "$bam" == 1 ]]; then
            dbam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_dedup.bam")
            sbam="$OUTPUT/$(basename "$(dirname "$(readlink -f "$dbam")")")/$(basename "${dbam%_*}")"
            samtools sort -@ "$THREADS" "$dbam" -o "${sbam}_sorted.bam"
            if [[ "$REMOVE" == true ]]; then
                rm "$dbam"
            fi
        elif [[ "$bam" == 0 ]]; then
            echo "Error: Filtered BAM files are required for the sorting step."
            exit 1
        fi
    done
    echo
    echo "Performing quality check."
    for i in "$OUTPUT"/*; do
        basename "$i"
        bam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam" | wc -l)
        if [[ "$bam" == 1 ]]; then
            sbam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam")
            sstats="$OUTPUT/$(basename "$(dirname "$(readlink -f "$sbam")")")/$(basename "${sbam%%.*}")"
            samtools stats -@ "$THREADS" "$sbam" > "${sstats}_stats.txt"
        fi
    done
    echo
    echo "Indexing reads."
    for i in "$OUTPUT"/*; do
        basename "$i"
        bam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam" | wc -l)
        if [[ "$bam" == 1 ]]; then
            sbam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam")
            sstats="$OUTPUT/$(basename "$(dirname "$(readlink -f "$sbam")")")/$sbam"
            samtools index -@ "$THREADS" -b "$sbam" "${sbam}.bai"
        fi
    done
    echo
}

function mapping() {
    if [[ "$STEP" == "mapping" ]]; then
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
        DIR="$INPUT"
    else
        DIR="$OUTPUT"
    fi
    echo "Calculating genome-wide coverage at each base pair."
    for i in "$DIR"/*; do
        basename "$i"
        bam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam" | wc -l)
        if [[ "$bam" == 1 ]]; then
            sbam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam")
            sbed="$OUTPUT/$(basename "$(dirname "$(readlink -f "$sbam")")")/$(basename "${sbam%%.*}")"
            samtools depth -a -o "${sbed}.bed" "$sbam"
        elif [[ "$bam" == 0 ]]; then
            echo "Error: Sorted BAM files are required for the mapping step."
            exit 1
        fi
    done
    echo
}
###############################################################################################################################################################
function peak_calling() {
    if [[ "$STEP" == "peak_calling" ]]; then
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
        DIR="$INPUT"
    else
        DIR="$OUTPUT"
    fi
    echo "Performing peak-calling."
    for i in "$DIR"/*"$ANTIBODY1"*; do
        basename "$i"
        sbam1=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam" -and -name "$ANTIBODY1" | wc -l)
        if [[ "$sbam1" == 1 ]]; then
            sbam1=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam" -and -name "$ANTIBODY1")
            npeaks="$(basename "${sbam1%_*}")"
            gsize=$(awk '{sum+=$2} END {print sum}' "$fa".fai)
            echo "$npeaks"
            echo "$gsize"
            exit1
            sbam2=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam" -and -name "$ANTIBODY2")
            macs2 callpeak -t "$sbam1" -c "$sbam2" -n "$npeaks"
        #elif [[ "$sbam" == 0 ]]; then
        #    echo "Error: Sorted BAM files are required for the mapping step."
        #    exit 1
        fi
    done
    echo
}


###################
## Run pipeline. ##
###################
#for s in "${FILTERED_STEPS[@]}"; do
#    "$s"
#done
peak_calling