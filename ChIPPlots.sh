#!/bin/bash

## Created: August 31, 2023
## Updated: August 31, 2023
## Author(s): Todd Lenz, tlenz001@ucr.edu
## ChIPPlots: Generates various plots for a ChIP-seq dataset.

function help {
    echo "ChIPPeaks.sh --help"
    echo "usage : ChIPPeaks.sh -m METADATA -c CONTROL -g GENOME [-h]"
    echo
    echo "---------------------------------------------------------------"
    echo " Input arguments:"
    echo "  -m|--metadata METADATA : ChIP-seq sample metadata file."
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
        "--genome") set -- "$@" "-g" ;;
        "--help") set -- "$@" "-h" ;;
        *) set -- "$@" "$arg"
    esac
done

pipedir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
SCRIPTS="$pipedir"/scripts

while getopts ":m:g:h:" opt; do
    case $opt in
        m) METADATA="$OPTARG";;
        g) GENOME="$OPTARG";;
        h) help;;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done

if [[ -z "$METADATA" ]]; then
    echo "Error: Script must be run with a metadata file."
    exit 1
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

for pkg in "pandas" "numpy" "matplotlib"; do
    if ! python3 -c "import $pkg" &> /dev/null; then
        echo "Installing $pkg."
        pip3 install "$pkg"
    fi
done

#######################################
## Check for necessary genome files. ##
#######################################
if [[ -n $GENOME ]]; then
    gff=$(find -L "$(readlink -f "$GENOME")" -mindepth 1 -maxdepth 1 \( -name "*.gff" \))
fi

#############################################
## Normalize samples and subtract control. ##
#############################################
data_dir="$(dirname "$(readlink -f "$METADATA")")"
if [[ "$(find -L "$(dirname "$(readlink -f "$METADATA")")" -mindepth 1 -maxdepth 3 -name "*_diff.bed" | wc -l)" -lt "$(("$(wc -l < "$METADATA")" - 1))" ]]; then
    echo "Normalizing samples and subtracting input/IGG."
    while IFS= read -r line; do
        SampleID="$(echo "$line" | awk -F '\t' '{print $1}')"
        echo "$SampleID"
        bed1="$(echo "$line" | awk -F '\t' '{print $5}')"
        bed1="$data_dir/${bed1#"$data_dir"}"
        bed1="${bed1%.*}.bed"
        bed2="$(echo "$line" | awk -F '\t' '{print $7}')"
        bed2="$data_dir/${bed2#"$data_dir"}"
        bed2="${bed2%.*}.bed"
        python3 "$SCRIPTS"/subtract_control.py -b1 "$bed1" -b2 "$bed2" -n CPM
    done < <(tail -n +2 "$METADATA")
    echo
fi

#############################################
## Normalize samples and subtract control. ##
#############################################
echo "Plotting genome-wide coverage."
while IFS= read -r line; do
    SampleID="$(echo "$line" | awk -F '\t' '{print $1}')"
    echo "$SampleID"
    bamReads="$(echo "$line" | awk -F '\t' '{print $5}')"
    bamReads="$data_dir/${bamReads#"$data_dir"}"
    bamReads="${bamReads%.*}_diff.bed"
    #python3 "$SCRIPTS"/plot_chromosome_coverage.py -b "$bamReads"
done < <(tail -n +2 "$METADATA")
echo


