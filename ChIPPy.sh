#!/bin/bash

## Created: June 13, 2022
## Updated: August 30, 2023
## Author(s): Todd Lenz, tlenz001@ucr.edu
## ChIPPy: A complete pipeline for ChIP-seq data analysis and plotting.


function help {
    echo "ChIPPy.sh --help"
    echo "usage : ChIPPy.sh -i INPUT -g GENOME [-o OUTPUT] [-s STEP] [-q QUALITY] [-t THREADS] [-r] [-h]"
    echo
    echo "----------------------------------------------------------------"
    echo " Required inputs:"
    echo "  -i|--input  INPUT        : Input data folder."
    echo "  -g|--genome GENOME       : Path to genome files."
    echo
    echo " Optional inputs:"
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
        "--treatment") set -- "$@" "-t" ;;
        "--control") set -- "$@" "-c" ;;
        "--threads") set -- "$@" "-p" ;;
        "--remove") set -- "$@" "-r" ;;
        "--help") set -- "$@" "-h" ;;
        *) set -- "$@" "$arg"
    esac
done

pipedir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
SCRIPTS="$pipedir"/scripts

while getopts ":i:g:o:s:q:t:r:h:" opt; do
    case $opt in
        i) INPUT="$OPTARG";;
        g) GENOME="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        s) STEP="$OPTARG";;
        q) QUALITY="$OPTARG";;
        t) THREADS="$OPTARG";;
        r) REMOVE=true;;
        h) help;;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done


############################
## Check input arguments. ##
############################
if [[ -z "$INPUT" ]]; then
    echo "Error: Input argument required."
    exit 1
elif [[ ! -e "$INPUT" ]]; then
    echo "Error: Input folder not found."
    exit 1
elif [[ -z "$GENOME" ]]; then
    echo "Error: Genome argument required."
    exit 1
elif [[ ! -e "$GENOME" ]]; then
    echo "Error: Genome folder not found."
    exit 1
elif [[ -n "$THREADS" && "$THREADS" -gt $(nproc) ]]; then
    echo "Error: Invalid value for threads argument."
    exit 1
elif [[ -n "$THREADS" && ! "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid value for threads argument."
    exit 1
fi


####################
## Define output. ##
####################
if [[ -e "$OUTPUT" && -z "$STEP" ]]; then
    echo "$OUTPUT folder alreads exists. Do you want to overwrite it? (y/n)"
    read -r ans
    if [[ "$ans" = "y" ]]; then
        rm -rf "$OUTPUT"
        mkdir "$OUTPUT"
        mkdir "$OUTPUT"/output_files
        OUTPUT="$OUTPUT"/output_files
        for i in "$INPUT"/*; do
            mkdir "$OUTPUT"/"$(basename "$i")"
        done
    fi
elif [[ -e "$OUTPUT" && -n "$STEP" ]]; then
    OUTPUT="$OUTPUT"/output_files
elif [[ -n "$OUTPUT" && -z "$STEP" ]]; then
    mkdir "$OUTPUT"
    mkdir "$OUTPUT"/output_files
    OUTPUT="$OUTPUT"/output_files
    for i in "$INPUT"/*; do
        mkdir "$OUTPUT"/"$(basename "$i")"
    done
elif [[ -n "$OUTPUT" && -n "$STEP" ]]; then
    OUTPUT="$OUTPUT"/output_files
elif [[ -z "$OUTPUT" ]]; then
    echo "No output folder selected. Files will be saved to input directory."
    mkdir "$INPUT"/output_files
    OUTPUT="$INPUT"/output_files
    for i in "$INPUT"/*; do
        mkdir "$OUTPUT"/"$(basename "$i")"
    done
fi


############################
## Define logs directory. ##
############################
if [[ ! -e "$(dirname "$OUTPUT")"/logs ]]; then
    log_dir="$(dirname "$OUTPUT")"/logs
    mkdir "$log_dir"
else
    log_dir="$(dirname "$OUTPUT")"/logs
fi


######################################################
## Check Python installation and required packages. ##
######################################################
for p in "python3" "fastqc" "bowtie2" "samtools"; do
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

if ! command -v pip3 &> /dev/null; then
    python -m ensurepip --upgrade
fi

if ! command -v picard &> /dev/null; then
    if [[ $(find -L "$SCRIPTS" -mindepth 1 -maxdepth 4 -name "picard.jar" | wc -l) -eq 0 ]]; then
        if [ ! -d "$SCRIPTS"/picard ]; then
            echo "Installing picard."
            git clone https://github.com/broadinstitute/picard.git "$SCRIPTS"/picard
        fi
        "$SCRIPTS"/picard/gradlew -p "$SCRIPTS"/picard shadowJar
    fi
fi


###############################################
## Define Phred score for quality filtering. ##
###############################################
if [[ -z $QUALITY ]]; then
    QUALITY=30
fi


#################################################
## Define number of processing threads to use. ##
#################################################
if [[ -z $THREADS ]]; then
    THREADS=$(($(nproc) / 2))
fi


###########################################
## Check for all necessary genome files. ##
###########################################
if [[ -n $GENOME ]]; then
    GENOME="$(readlink -f "$GENOME")"
    fa=$(find -L "$GENOME" -mindepth 1 -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \))
    bname="${fa%%.*}"
    if [[ $fa = "" && ! -e "$bname.1.bt2" ]]; then
        echo "Error: No FASTA found to build Bowtie2 indexes."
    elif [[ ! -e "$bname.1.bt2" ]]; then
        echo "Creating Bowtie2 indexes."
        bowtie2-build --threads "$THREADS" "$fa" "$bname"
    else
        echo "Bowtie2 indexes detected."
    fi
    if [[ $fa = "" && ! -e "$fa.fai" ]]; then
        echo "Error: No FASTA found to build FASTA index."
    elif [[ ! -e "$fa.fai" ]]; then
        echo "Creating FASTA index."
        samtools faidx "$fa"
    else
        echo "FASTA index detected."
    fi
fi


####################################################
## Check starting step and required files option. ##
####################################################
STEP="${STEP/^ //}"
ALL_STEPS=("quality_check" "trimming" "alignment" "deduplication" "filtering" "sorting" "mapping")
QUALITY_CHECK="quality_check"
TRIMMING="trimming"
ALIGNMENT="alignment"
DEDUPLICATION="deduplication"
FILTERING="filtering"
SORTING="sorting"
MAPPING="mapping"

if [[ -n $STEP ]]; then
    case $STEP in
        "quality_check") FILTERED_STEPS=("${ALL_STEPS[@]}");;
        "trimming") FILTERED_STEPS=("${ALL_STEPS[@]:1}");;
        "alignment") FILTERED_STEPS=("${ALL_STEPS[@]:2}");;
        "deduplication") FILTERED_STEPS=("${ALL_STEPS[@]:3}");;
        "filtering") FILTERED_STEPS=("${ALL_STEPS[@]:4}");;
        "sorting") FILTERED_STEPS=("${ALL_STEPS[@]:5}");;
        "mapping") FILTERED_STEPS=("${ALL_STEPS[@]:6}");;
        *) echo "Invalid step: $STEP"; exit 1;;
    esac
else
    FILTERED_STEPS=("${ALL_STEPS[@]}")
fi

if [[ -n "$STEP" ]]; then
    for s in $STEP; do
        for i in "${QUALITY_CHECK[@]}"; do
            if [[ "$i" = "$s" ]]; then
                fq=$(find -L "$INPUT" -mindepth 2 -maxdepth 2  -name "*.fastq*" -o -name "*.fastq*" | wc -l)
                if [[ "$fq" -eq 0 ]]; then
                    echo "Error: No '.fastq(.gz)' files detected."
                    exit 1
                elif [[ "$fq" -eq 2 ]]; then
                    fq_r1=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R1*" | wc -l)
                    fq_r2=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R2*" | wc -l)
                    if [[ "$fq_r1" -eq 0 || "$fq_r2" -eq 0 || "$fq_r1" -ne "$fq_r2" ]]; then
                        echo "Error: Paired-end input files must contain '_R1' and '_R2' in the file names."
                        exit 1
                    fi
                fi
            fi
        done
        for i in "${TRIMMING[@]}"; do
            if [[ "$i" = "$s" ]]; then
                fq=$(find -L "$INPUT" -mindepth 2 -maxdepth 2  -name "*.fastq*" -o -name "*.fastq*" | wc -l)
                if [[ "$fq" -eq 0 ]]; then
                    echo "Error: No '.fastq(.gz)' files detected."
                    exit 1
                elif [[ "$fq" -eq 2 ]]; then
                    fq_r1=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R1*" | wc -l)
                    fq_r2=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R2*" | wc -l)
                    if [[ "$fq_r1" -eq 0 || "$fq_r2" -eq 0 || "$fq_r1" -ne "$fq_r2" ]]; then
                        echo "Error: Paired-end input files must contain '_R1' and '_R2' in the file names."
                        exit 1
                    fi
                fi
            fi
        done
        for i in "${ALIGNMENT[@]}"; do
            if [[ "$i" = "$s" ]]; then
                tfq=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 -name "*_paired.fastq*" -o -name "*_paired.fq*" -o -name "*_trimmed.fastq*" -o -name "*_trimmed.fq*" | wc -l)
                if [[ "$tfq" -eq 0 ]]; then
                    echo "Error: Trimmed fastq files are required to start the pipeline at the alignment step."
                    exit 1
                elif [[ "$tfq" -eq 1 ]]; then
                    tfq_r1=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*_trimmed.fastq*" -o -name "*_trimmed.fq*" \) -and -name "*_R1*" | wc -l)
                    if [[ "$tfq_r1" -eq 0 ]]; then
                        echo "Error: Trimmed fastq files are required to start the pipeline at the alignment step."
                        echo "Ensure that single-end input files contain '_trimmed.fastq' and '_R1' in the file name."
                        exit 1
                    fi
                elif [[ "$tfq" -eq 2 ]]; then
                    tfq_r1=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*_paired.fastq*" -o -name "*_paired.fq*" \) -and -name "*_R1*" | wc -l)
                    tfq_r2=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 \( -name "*_paired.fastq*" -o -name "*_paired.fq*" \) -and -name "*_R2*" | wc -l)
                    if [[ "$tfq_r1" -eq 0 || "$tfq_r2" -eq 0 || "$tfq_r1" -ne "$tfq_r2" ]]; then
                        echo "Error: Trimmed fastq files are required to start the pipeline at the alignment step."
                        echo "Ensure that paired-end input files contain '_paired.fastq' and '_R1'/'_R2' in the file name."
                        exit 1
                    fi
                fi
            fi
        done
        for i in "${DEDUPLICATION[@]}"; do
            if [[ "$i" = "$s" ]]; then
                asam=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 -name "*_aligned.sam" | wc -l)
                if [[ "$asam" -eq 0 ]]; then
                    echo "Error: Aligned SAM files are required to start the pipeline at the deduplication step."
                    exit 1
                fi
            fi
        done
        for i in "${FILTERING[@]}"; do
            if [[ "$i" = "$s" ]]; then 
                dsam=$(find -L "$INPUT" -mindepth 2 -maxdepth 2 -name "*_dedup.sam" | wc -l)
                if [[ "$dsam" -eq 0 ]]; then
                    echo "Error: Deduplicated SAM files are required to start the pipeline at the filtering step."
                    exit 1
                fi
            fi
        done
        for i in "${SORTING[@]}"; do
            if [[ "$i" = "$s" ]]; then 
                fbam=$(find -L "$INPUT" -mindepth 2 -maxdepth 4 -name "*_dedup.bam" | wc -l)
                if [[ "$fbam" -eq 0 ]]; then
                    echo "Error: Filtered BAM files are required to start the pipeline at the sorting step."
                    exit 1
                fi
            fi
        done
        for i in "${MAPPING[@]}"; do
            if [[ "$i" = "$s" ]]; then 
                sbam=$(find -L "$INPUT" -mindepth 2 -maxdepth 4 -name "*_sorted.bam" | wc -l)
                if [[ "$sbam" -eq 0 ]]; then
                    echo "Error: Sorted BAM files are required to start the pipeline at the mapping step."
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


###########################################
## Functions/steps in ChIP-seq pipeline. ##
###########################################
function quality_check() {
    echo "Performing FastQC quality check."
    echo "----------------------------------------------------------------"
    for i in "$INPUT"/*; do
        basename "$i"
        fq=$(find -L "$i" -mindepth 1 -maxdepth 1  -name "*.fastq*" -o -name "*.fastq*" | wc -l)
        if [[ "$fq" -eq 1 ]]; then
            fq_r1=$(find -L "$i" -mindepth 1 -maxdepth 1  -name "*.fastq*" -o -name "*.fastq*" )
            fastqc -t "$THREADS" "$fq_r1"
            if [[ "$OUTPUT" != "$INPUT" ]]; then
                mv "${fq_r1%%.*}"_fastqc.html "$OUTPUT/$(basename "$(dirname "$fq_r1")")/$(basename "${fq_r1%%.*}")_fastqc.html"
            fi
            rm "${fq_r1%%.*}"_fastqc.zip
        fi
        if [[ "$fq" -eq 2 ]]; then
            fq_r1=$(find -L "$i" -mindepth 1 -maxdepth 1 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R1*")
            fq_r2=$(find -L "$i" -mindepth 1 -maxdepth 1 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R2*")
            fastqc -t "$THREADS" "$fq_r1" >> "$log_dir"/fastqc.log 2>&1
            fastqc -t "$THREADS" "$fq_r2" >> "$log_dir"/fastqc.log 2>&1
            if [[ "$OUTPUT" != "$INPUT" ]]; then
                mv "${fq_r1%%.*}"_fastqc.html "$OUTPUT/$(basename "$(dirname "$fq_r1")")/$(basename "${fq_r1%%.*}")_fastqc.html"
                mv "${fq_r2%%.*}"_fastqc.html "$OUTPUT/$(basename "$(dirname "$fq_r2")")/$(basename "${fq_r2%%.*}")_fastqc.html"
            fi
            rm "${fq_r1%%.*}"_fastqc.zip
            rm "${fq_r2%%.*}"_fastqc.zip
        elif [[ "$fq" -eq 0 ]]; then
            echo "Error: Input files must be '.fastq(.gz)'."
            exit 1
        fi
    done
    echo
}

function trimming() {
    if [[ "$STEP" = "trimming" ]]; then
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
    fi
    echo "Trimming and pairing R1 and R2 reads."
    echo "----------------------------------------------------------------"
    for i in "$INPUT"/*; do
        basename "$i"
        fq=$(find -L "$i" -mindepth 1 -maxdepth 1  -name "*.fastq*" -o -name "*.fastq*" | wc -l)
        if ! command -v trimmomatic &> /dev/null; then
            tmjar=$(find -L "$SCRIPTS" -mindepth 1 -maxdepth 3 -name "*trimmomatic*")
        fi
        if [[ "$fq" -eq 1 ]]; then
            fq_r1=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*.fastq*" -o -name "*.fq*" )
            out_r1="$OUTPUT"/"$(basename "$(dirname "$fq_r1")")"/"$(basename "${fq_r1%%.*}")"
            if ! command -v trimmomatic &> /dev/null; then
                java -jar "$tmjar" SE -threads "$THREADS" "$fq_r1" "${out_r1}_trimmed.fastq.gz" LEADING:"$QUALITY" TRAILING:"$QUALITY" MINLEN:25 ILLUMINACLIP:"$GENOME"/adapters.txt:2:30:10 >> "$log_dir"/trimming.log 2>&1
            else
                trimmomatic SE -threads "$THREADS" "$fq_r1" "${out_r1}_trimmed.fastq.gz" LEADING:"$QUALITY" TRAILING:"$QUALITY" MINLEN:25 ILLUMINACLIP:"$GENOME"/adapters.txt:2:30:10 >> "$log_dir"/trimming.log 2>&1
            fi
        elif [[ "$fq" -eq 2 ]]; then
            fq_r1=$(find -L "$i" -mindepth 1 -maxdepth 1 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R1*")
            fq_r2=$(find -L "$i" -mindepth 1 -maxdepth 1 \( -name "*.fastq*" -o -name "*.fq*" \) -and -name "*_R2*")
            out_r1="$OUTPUT"/"$(basename "$(dirname "$fq_r1")")"/"$(basename "${fq_r1%%.*}")"
            out_r2="$OUTPUT"/"$(basename "$(dirname "$fq_r2")")"/"$(basename "${fq_r2%%.*}")"
            if ! command -v trimmomatic &> /dev/null; then
                java -jar "$tmjar" PE -threads "$THREADS" "$fq_r1" "$fq_r2" "${out_r1}_paired.fastq.gz" "${out_r1}_unpaired.fastq.gz" "${out_r2}_paired.fastq.gz" "${out_r2}_unpaired.fastq.gz" ILLUMINACLIP:"$GENOME"/adapters.txt:2:30:10 LEADING:"$QUALITY" TRAILING:"$QUALITY" MINLEN:25 >> "$log_dir"/trimming.log 2>&1
            else
                trimmomatic PE -threads "$THREADS" "$fq_r1" "$fq_r2" "${out_r1}_paired.fastq.gz" "${out_r1}_unpaired.fastq.gz" "${out_r2}_paired.fastq.gz" "${out_r2}_unpaired.fastq.gz" ILLUMINACLIP:"$GENOME"/adapters.txt:2:30:10 LEADING:"$QUALITY" TRAILING:"$QUALITY" MINLEN:25 >> "$log_dir"/trimming.log 2>&1
            fi
        elif [[ "$fq" -eq 0 ]]; then
            echo "Error: Input files must be '.fastq(.gz)'."
            exit 1
        fi
    done
    echo
}

function alignment() {
    if [[ "$STEP" = "alignment" ]]; then
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
    fi
    if [[ ! -e "$bname.1.bt2" ]]; then
        echo "Error: No Bowtie2 indexes found."
        exit 1
    fi
    echo "Aligning reads to the genome."
    echo "----------------------------------------------------------------"
    for i in "$OUTPUT"/*; do
        basename "$i"
        tfq=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_trimmed.fastq*" -o -name "*_trimmed.fq*" -o -name "*_paired.fastq*" -o -name "*_paired.fq*"  | wc -l)
        fa="$(find -L "$(readlink -f "$GENOME")" -mindepth 1 -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \))"
        bidx="${fa%%.*}"
        if [[ $tfq -eq 1 ]]; then
            tfq_r1=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_trimmed.fastq*" -o -name "*_trimmed.fq*" )
            out="$OUTPUT/$(basename "$(dirname "$(readlink -f "$tfq_r1")")")/$(basename "${tfq_r1%_*}")"
            bowtie2 -x "$bidx" -1 "$tfq_r1" -S "${out}_aligned.sam" -p "$THREADS" --very-sensitive >> "$log_dir"/alignment.log 2>&1
            if [[ "$REMOVE" = true ]]; then
                rm "$i"/*_trimmed.fastq*
            fi
        elif [[ $tfq -eq 2 ]]; then
            tfq_r1=$(find -L "$i" -mindepth 1 -maxdepth 1 \( -name "*_paired.fastq*" -o -name "*_paired.fq*" \) -and -name "*_R1*")
            tfq_r2=$(find -L "$i" -mindepth 1 -maxdepth 1 \( -name "*_paired.fastq*" -o -name "*_paired.fq*" \) -and -name "*_R2*")
            out="$OUTPUT/$(basename "$(dirname "$(readlink -f "$tfq_r1")")")/$(basename "${tfq_r1%_*}")"
            asam="${out/_R1/}"
            bowtie2 -x "$bidx" -1 "$tfq_r1" -2 "$tfq_r2" -S "${asam}_aligned.sam" -p "$THREADS" --very-sensitive >> "$log_dir"/alignment.log 2>&1
            if [[ "$REMOVE" = true ]]; then
                rm "$i"/*_paired.fastq* "$i"/*_unpaired.fastq*
            fi
        elif [[ "$tfq" -eq 0 ]]; then
            echo "Error: Trimmed '.fastq(.gz)' files are required for the alignment step."
            exit 1
        fi
    done
    echo
}

function deduplication() {
    if [[ "$STEP" = "deduplication" ]]; then
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
    fi
    echo "Removing PCR duplicates."
    echo "----------------------------------------------------------------"
    for i in "$OUTPUT"/*; do
        basename "$i"
        sam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_aligned.sam" | wc -l)
        if ! command -v picard &> /dev/null; then
            ptjar=$(find -L "$SCRIPTS" -mindepth 1 -maxdepth 4 -name "picard.jar")
        fi
        if [[ "$sam" -eq 1 ]]; then
            asam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_aligned.sam")
            dsam="$OUTPUT/$(basename "$(dirname "$(readlink -f "$asam")")")/$(basename "${asam%_*}")"
            if ! command -v picard &> /dev/null; then
                java -jar "$ptjar" MarkDuplicates -I "$asam" -O "${dsam}_dedup.sam" -M "${dsam}_dedup_metrics.txt" -ASO queryname >> "$log_dir"/deduplication.log 2>&1
            else
                picard MarkDuplicates -I "$asam" -O "${dsam}_dedup.sam" -M "${dsam}_dedup_metrics.txt" -ASO queryname >> "$log_dir"/deduplication.log 2>&1
            fi
            if [[ "$REMOVE" = true ]]; then
                rm "$asam"
            fi
        elif [[ "$sam" -eq 0 ]]; then
            echo "Error: SAM files are required for the deduplication step."
            exit 1
        fi
    done
    echo
}

function filtering() {
    if [[ "$STEP" = "filtering" ]]; then
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
    fi
    echo "Filtering low quality reads."
    echo "----------------------------------------------------------------"
    for i in "$OUTPUT"/*; do
        basename "$i"
        sam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_dedup.sam" | wc -l)
        if [[ "$sam" -eq 1 ]]; then
            dsam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_dedup.sam")
            dbam="$OUTPUT/$(basename "$(dirname "$(readlink -f "$dsam")")")/$(basename "${dsam%%.*}")"
            samtools view -@ "$THREADS" -q "$QUALITY" -f 0X02 -F 0X04 -b "$dsam" -o "${dbam}.bam" >> "$log_dir"/filtering.log 2>&1
            if [[ "$REMOVE" = true ]]; then
                rm "$dsam"
            fi
        elif [[ "$sam" -eq 0 ]]; then
            echo "Error: Deduplicated SAM files are required for the filtering step."
            exit 1
        fi
    done
    echo
}

function sorting() {
    if [[ "$STEP" = "sorting" ]]; then
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
    fi
    echo "Sorting reads."
    echo "----------------------------------------------------------------"
    for i in "$OUTPUT"/*; do
        basename "$i"
        bam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_dedup.bam" | wc -l)
        if [[ "$bam" -eq 1 ]]; then
            dbam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_dedup.bam")
            sbam="$OUTPUT/$(basename "$(dirname "$(readlink -f "$dbam")")")/$(basename "${dbam%_*}")"
            samtools sort -@ "$THREADS" "$dbam" -o "${sbam}_sorted.bam" >> "$log_dir"/sorting.log 2>&1
            if [[ "$REMOVE" = true ]]; then
                rm "$dbam"
            fi
        elif [[ "$bam" -eq 0 ]]; then
            echo "Error: Filtered BAM files are required for the sorting step."
            exit 1
        fi
    done
    echo
    echo "Performing quality check."
    echo "----------------------------------------------------------------"
    for i in "$OUTPUT"/*; do
        basename "$i"
        bam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam" | wc -l)
        if [[ "$bam" -eq 1 ]]; then
            sbam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam")
            sstats="$OUTPUT/$(basename "$(dirname "$(readlink -f "$sbam")")")/$(basename "${sbam%%.*}")"
            samtools stats -@ "$THREADS" -d -x "$sbam" > "${sstats}_stats.txt"
        fi
    done
    echo
    echo "Indexing reads."
    echo "----------------------------------------------------------------"
    for i in "$OUTPUT"/*; do
        basename "$i"
        bam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam" | wc -l)
        if [[ "$bam" -eq 1 ]]; then
            sbam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam")
            sstats="$OUTPUT/$(basename "$(dirname "$(readlink -f "$sbam")")")/$sbam"
            samtools index -@ "$THREADS" -b "$sbam" "${sbam}.bai"
        fi
    done
    echo
}

function mapping() {
    if [[ "$STEP" = "mapping" ]]; then
        echo "Starting ChIP-seq analysis pipeline at the $STEP step."
        echo
    fi
    echo "Calculating genome-wide coverage at each base pair."
    echo "----------------------------------------------------------------"
    for i in "$OUTPUT"/*; do
        basename "$i"
        bam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam" | wc -l)
        if [[ "$bam" -eq 1 ]]; then
            sbam=$(find -L "$i" -mindepth 1 -maxdepth 1 -name "*_sorted.bam")
            sbed="$OUTPUT/$(basename "$(dirname "$(readlink -f "$sbam")")")/$(basename "${sbam%%.*}")"
            samtools depth -a -o "${sbed}.bed" "$sbam" >> "$log_dir"/mapping.log 2>&1
        elif [[ "$bam" -eq 0 ]]; then
            echo "Error: Sorted BAM files are required for the mapping step."
            exit 1
        fi
    done
    echo
}


###################
## Run pipeline. ##
###################
for fs in "${FILTERED_STEPS[@]}"; do
    "$fs"
done
