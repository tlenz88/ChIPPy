#!/bin/bash

## Created: August 25, 2023
## Updated: August 28, 2023
## Author(s): Todd Lenz, tlenz001@ucr.edu
## ChIPPeaks: Performs peak calling and differential peak calling using an input chip-seq metadata file describing samples. See example metadata file.


function help {
    echo "ChIPPeaks.sh --help"
    echo "usage : ChIPPeaks.sh [-g GENOME] [-m METADATA] [-o OUTPUT] [-h]"
    echo
    echo "-----------------------------------------------------------------------"
    echo " Optional inputs:"
    echo "  -g|--genome GENOME     : Path to genome files."
    echo "  -m|--metadata METADATA : ChIP-seq sample metadata file."
    echo "  -o|--output OUTPUT     : Output directory."
    echo "  -h|--help HELP         : Show help message."
    echo "-----------------------------------------------------------------------"
    exit 0;
}


#############################
## Define input arguments. ##
#############################
for arg in "$@"; do
    case "$arg" in
        "--genome") set -- "$@" "-g" ;;
        "--metadata") set -- "$@" "-m" ;;
        "--output") set -- "$@" "-o" ;;
        "--help") set -- "$@" "-h" ;;
        *) set -- "$@" "$arg"
    esac
done

pipedir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
SCRIPTS="$pipedir"/scripts

while getopts ":g:m:o:h:" opt; do
    case $opt in
        g) GENOME="$OPTARG";;
        m) METADATA="$OPTARG";;
        o) OUTPUT="$OPTARG";;
        h) help;;
        *) echo "Error: '$OPTARG' is an invalid argument."
    esac
done


##########################################
## Define input and output directories. ##
##########################################
if [[ -z "$OUTPUT" && $METADATA != "" ]]; then
    OUTPUT="$(awk 'NR==2 {split($5, path_parts, "/"); path = ""; for (i=1; i<length(path_parts)-2; i++) path = path path_parts[i] "/"; print path}' $METADATA)"
elif [[ -z "$OUTPUT" && -z "$METADATA" ]]; then
    OUTPUT="$pipedir"
fi
OUTPUT="$(basename $OUTPUT)"

if [[ -z "$GENOME" ]]; then
    GENOME="$pipedir"/genomes
fi

if [[ -z "$METADATA" ]]; then
    METADATA="$pipedir"/examples/example_metadata.txt
fi

if [[ ! -e "$OUTPUT"/figures ]]; then
    fig_dir="$OUTPUT"/figures
    mkdir "$fig_dir"
else
    fig_dir="$OUTPUT"/figures
fi

if [[ ! -e "$OUTPUT"/logs ]]; then
    log_dir="$OUTPUT"/logs
    mkdir "$log_dir"
else
    log_dir="$OUTPUT"/logs
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


##############################
## Run peak calling script. ##
##############################
python3 "$SCRIPTS"/peak_calling.py -g $GENOME -m $METADATA -f $fig_dir -l $log_dir
