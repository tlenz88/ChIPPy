#!/bin/bash

## Created: August 25, 2023
## Updated: August 31, 2023
## Author(s): Todd Lenz, tlenz001@ucr.edu
## ChIPPeaks: Performs peak calling and differential peak calling using an input chip-seq metadata file describing samples.


function help {
    echo "ChIPPeaks.sh --help"
    echo "usage : ChIPPeaks.sh -m METADATA -c CONTROL -g GENOME [-h]"
    echo
    echo "---------------------------------------------------------------"
    echo " Input arguments:"
    echo "  -m|--metadata METADATA : ChIP-seq sample metadata file."
    echo "  -c|--control CONTROL   : Control condition."
    echo "  -g|--genome GENOME     : Directory containing genome files."
    echo "  -h|--help HELP         : Show help message."
    echo "---------------------------------------------------------------"
    exit 0;
}


#############################
## Define input arguments. ##
#############################
for arg in "$@"; do
    case "$arg" in
        "--metadata") set -- "$@" "-m" ;;
        "--control") set -- "$@" "-c" ;;
        "--genome") set -- "$@" "-g" ;;
        "--help") set -- "$@" "-h" ;;
        *) set -- "$@" "$arg"
    esac
done

pipedir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
SCRIPTS="$pipedir"/scripts

while getopts ":m:c:g:h:" opt; do
    case $opt in
        m) METADATA="$OPTARG";;
        c) CONTROL="$OPTARG";;
        g) GENOME="$OPTARG";;
        h) help;;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done

if [[ -z "$METADATA" ]]; then
    echo "Error: Script must be run with a metadata file."
    exit 1
fi


########################################
## Define figure and log directories. ##
########################################
data_dir="$(dirname "$(readlink -f "$METADATA")")"

if [[ ! -e "$data_dir"/figures ]]; then
    fig_dir="$data_dir"/figures
    mkdir "$fig_dir"
else
    fig_dir="$data_dir"/figures
fi

if [[ ! -e "$data_dir"/logs ]]; then
    log_dir="$data_dir"/logs
    mkdir "$log_dir"
else
    log_dir="$data_dir"/logs
fi

######################################################
## Check Python installation and required packages. ##
######################################################
if ! command -v "python3" &> /dev/null; then
    echo "Installing newest version of python3."
    if command -v apt &> /dev/null; then
        sudo apt update
        sudo apt upgrade
        sudo apt install -y python3
    elif command -v yum &> /dev/null; then
        sudo yum install -y python3
    elif command -v dnf &> /dev/null; then
        sudo dnf update
        sudo dnf install -y python3
    elif command -v zypper &> /dev/null; then
        sudo zypper refresh
        sudo zypper update
        sudo zypper --non-interactive install python3
    elif command -v pacman &> /dev/null; then
        sudo pacman -Syu
        sudo pacman install --noconfirm python3
    elif command -v brew &> /dev/null; then
        brew update
        brew upgrade
        brew install -y python3
    else
        echo "Can't determine package manager. Install python3 manually."
        exit 1
    fi
fi

if ! command -v pip3 &> /dev/null; then
    python -m ensurepip --upgrade
fi

if ! command -v macs2 &> /dev/null; then
    echo "Installing macs2."
    pip3 install macs2
fi


#######################################
## Check for necessary genome files. ##
#######################################
function get_chrom_sizes() {
    fa=$(find -L "$1" -mindepth 1 -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \))
    chrom_sizes=$(find -L "$1" -mindepth 1 -maxdepth 1 \( -name "*.chrom.sizes" \))
    bname="${fa%%.*}"
    if [[ -e "$bname.chrom.sizes"  ]]; then
        echo "'.chrom.sizes' file found."
        gsize="$(awk -F '\t' '{ sum += $2 } END { print sum }' "$chrom_sizes")"
    elif [[ $fa = "" && ! -e "$bname.chrom.sizes"  ]]; then
        echo "Error: No FASTA or '.chrom.sizes' file found."
        exit 1
    elif [[ $fa != "" && ! -e "$bname.chrom.sizes" ]]; then
        echo "Finding chromosome lengths."
        python3 "$SCRIPTS"/create_sizes.py "$fa"
        gsize="$(awk -F '\t' '{ sum += $2 } END { print sum }' "$chrom_sizes")"
    fi
    echo
}

if [[ -z $GENOME ]]; then
    get_chrom_sizes "$pipedir"/genomes
else
    genome_dir="$(readlink -f "$GENOME")"
    get_chrom_sizes "$genome_dir"
fi


#######################
## Run peak calling. ##
#######################
function check_bam_type() {
    if $(samtools view -f 0x1 "$bamControl" | head -n 1 | wc -l) -eq 1; then
        echo "BAMPE"
    else
        echo "BAM"
    fi
}

function check_alg_type() {
    Factor=$(echo "$1" | awk -F '\t' '{print $2}')
    if [[ $Factor =~ ^H2[AB]K || $Factor =~ ^H3K ]]; then
        echo "--broad"
    else
        echo ""
    fi
}

echo "Performing peak calling."
while IFS= read -r line; do
    SampleID="$(echo "$line" | awk -F '\t' '{print $1}')"
    echo "$SampleID"
    bamReads="$(echo "$line" | awk -F '\t' '{print $5}')"
    bamReads="$data_dir/${bamReads#"$data_dir"}"
    bamControl="$(echo "$line" | awk -F '\t' '{print $7}')"
    bamControl="$data_dir/${bamControl#"$data_dir"}"
    bam_format="$(check_bam_type "$bamReads")"
    outdir="$(dirname "$bamReads")"
    alg="$(check_alg_type "$line")"
    macs2 callpeak -t "$bamReads" -c "$bamControl" -f "$bam_format" -g "$gsize" -n "$SampleID" -q 0.05 --outdir "$outdir" -B "$alg" >> "$log_dir"/peak_calling.log 2>&1
done < <(tail -n +2 "$METADATA")
echo

while IFS= read -r line; do
    if [[ "$(echo "$line" | awk -F'\t' '{print NF}')" -eq 7 ]]; then
        if [[ "$(echo "$line" | awk -F '\t' '{print $1}')" = "SampleID" ]]; then
            echo -e "$line\tPeaks\tPeakCaller"
        else
            bamReads="$(echo "$line" | awk -F '\t' '{print $5}')"
            echo -e "$line\t${bamReads%_*}_peaks.xls\tmacs"
        fi
    else
        echo "$line"
    fi
done < "$METADATA" > "${METADATA%.*}_updated.txt"
mv "${METADATA%.*}_updated.txt" "$METADATA"


####################################
## Run differential peak calling. ##
####################################
if [[ -z "$CONTROL" ]]; then
    echo "Warning: No control argument entered. Differential peak calling step can not be performed."
elif [[ $(awk '{ print $3 }' < "$METADATA" | grep -c "$CONTROL") -eq 0 ]]; then
    echo "Warning: Control argument not found in metadata. Differential peak calling step can not be performed."
else
    echo "Finding differential peaks."
    Rscript "$SCRIPTS"/differential_peak_calling.R "$(dirname "$METADATA")" "$(basename "$METADATA")" "$CONTROL" >> "$log_dir"/diff_peak_calling.log 2>&1
fi
