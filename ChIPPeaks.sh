#!/bin/bash

## Created: August 25, 2023
## Updated: August 28, 2023
## Author(s): Todd Lenz, tlenz001@ucr.edu
## ChIPPeaks: Performs peak calling and differential peak calling using an input chip-seq metadata file describing samples.


function help {
    echo "ChIPPeaks.sh --help"
    echo "usage : ChIPPeaks.sh -m METADATA -c CONTROL -g GENOME [-h]"
    echo
    echo "---------------------------------------------------------------"
    echo " Input arguments:"
    echo "  -m|--metadata METADATA : ChIP-seq sample metadata file."
    echo "  -c|--control CONTROL   : Control \'Condition\'."
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

while getopts ":m:c:o:g:h:" opt; do
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
data_dir="$(dirname "$(readlink -f $METADATA)")"

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
for p in "python3"; do
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

if ! command -v macs2 &> /dev/null; then
    echo "Installing macs2."
    pip3 install macs2
fi

for pkg in "pandas"; do
    if ! python3 -c "import $pkg" &> /dev/null; then
        echo "Installing $pkg."
        pip3 install "$pkg"
    fi
done


#######################################
## Check for necessary genome files. ##
#######################################
function get_chrom_sizes() {
    fa=$(find -L "$1" -mindepth 1 -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \))
    chrom_sizes=$(find -L "$1" -mindepth 1 -maxdepth 1 \( -name "*.chrom.sizes" \))
    bname="${fa%%.*}"
    if [[ $fa == "" && ! -e "$bname.chrom.sizes"  ]]; then
        echo "Error: No FASTA or '.chrom.sizes' file found."
    elif [[ $fa != "" && ! -e "$bname.chrom.sizes" ]]; then
        echo "Finding chromosome lengths."
        python3 "$SCRIPTS"/create_sizes.py "$fa"
    fi
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
    if $(samtools view -f 0x1 "$bamControl" | head -n 1 | wc -l) == 1; then
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

gsize="$(awk -F '\t' '{ sum += $2 } END { print sum }' $chrom_sizes)"

echo "Performing peak calling."
while IFS= read -r line; do
    bamReads="$(echo "$line" | awk -F '\t' '{print $5}')"
    bamReads="$data_dir/${bamReads#$data_dir}"
    bamControl="$(echo "$line" | awk -F '\t' '{print $7}')"
    bamControl="$data_dir/${bamControl#$data_dir}"
    bam_format="$(check_bam_type "$bamReads")"
    name="$(echo "$line" | awk -F '\t' '{print $1}')"
    echo "$name"
    outdir="$(dirname "$bamReads")"
    alg="$(check_alg_type "$line")"
    macs2 callpeak -t "$bamReads" -c "$bamControl" -f "$bam_format" -g "$gsize" -n "$name" -q 0.05 --outdir "$outdir" -B "$alg" >> "$log_dir"/peak_calling.log 2>&1
done < <(tail -n +2 "$METADATA")


####################################
## Run differential peak calling. ##
####################################
if [[ -z "$CONTROL" ]]; then
    echo "No control argument entered. Differential peak calling step can not be performed."
elif [[ $(awk '{ print $3 }' < "$METADATA" | grep "$CONTROL" | wc -l) == 0 ]]; then
    echo "Control argument not found metadata."
else
    echo "Finding differential peaks."
    Rscript "$SCRIPTS"/differential_peak_calling.R "$METADATA" "$CONTROL" >> "$log_dir"/diff_peak_calling.log 2>&1
fi
