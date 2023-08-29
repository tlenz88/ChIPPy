#!/usr/bin/env python3

"""
Created: August 23, 2023
Updated: August 28, 2023
Author(s): Todd Lenz, tlenz001@ucr.edu

Runs MACS2 callpeak using a ChIP-seq metadata file. An example ChIP-seq
metadata file is available in the 'examples' folder.
"""

import sys
import argparse
import os
import subprocess
import re
import glob
import pandas as pd
from Bio import SeqIO
import pysam


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-m',
                        '--metadata',
                        dest='metadata',
                        help='Path to ChIP-seq metadata file.',
                        required=True)
    parser.add_argument('-g',
                        '--genome',
                        dest='genome',
                        help='Path to directory containing genome files.',
                        required=False)
    parser.add_argument('-o',
                        '--out_dir',
                        dest='out_dir',
                        help='Path to output directory.',
                        required=False)
    parser.add_argument('-f',
                        '--fig_dir',
                        dest='fig_dir',
                        help='Path to figures directory.',
                        required=False)
    parser.add_argument('-l',
                        '--log_dir',
                        dest='log_dir',
                        help='Path to logs directory..',
                        required=False)
    return parser.parse_args()


def input_params(args):
    metadata = pd.read_csv(args.metadata, sep = '\t')
    try:
        genome_dir = os.path.abspath(args.genome)
    except:
        genome_dir = ''.join([os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), '/genomes/'])
    try:
        fig_dir = os.path.abspath(args.fig_dir)
    except:
        fig_dir = ''.join([os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), '/figures/'])
    try:
        log_dir = os.path.abspath(args.log_dir)
    except:
        log_dir = ''.join([os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), '/logs/'])
    return metadata, genome_dir, fig_dir, log_dir


def genome_files(genome_dir):
    cs = glob.glob(os.path.join(genome_dir, '*.chrom.sizes'))
    if len(cs) == 0:
        fa = glob.glob(os.path.join(genome_dir, '*.fasta'))
        if len(fa) == 1:
            for f in fa:
                fasta = f
        elif len(fa) == 0:
            print('Error: No FASTA or \'.chrom.sizes\' file found.')
            sys.exit()
        else:
            print('Error: Multiple FASTA files found \
                  in the \'genomes\' directory.')
            sys.exit()
        create_sizes = ''.join([os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), '/scripts/create_sizes.py'])
        subprocess.run(['python3', create_sizes, fasta])
        for s in glob.glob(os.path.join(genome_dir, '*.chrom.sizes')):
            return pd.read_csv(s, sep='\t', header=None)
    elif len(cs) == 1:
        for s in cs:
            print('\'.chrom.sizes\' file detected.')
            return pd.read_csv(s, sep='\t', header=None)
    else:
        print('Error: Multiple \'.chrom.sizes\' files \
              in \'genomes\' directory.')
        sys.exit()


def peak_calling(sample, gsize, log_dir):
    bamReads = os.path.abspath(sample[5])
    bamControl = os.path.abspath(sample[7])
    bam_format = check_bam_type(bamReads)
    name = str(sample[1])
    outdir = os.path.dirname(bamReads)
    alg = check_alg_type(sample[2])
    print(sample[1])
    with open(os.path.join(log_dir, 'peak_calling.log'), 'a+') as log_file:
        subprocess.run(['macs2', 'callpeak', '-t', bamReads, '-c', 
                        bamControl, '-f', bam_format, '-g', gsize, '-n', 
                        name, '-q', '0.05', '--outdir', outdir, '-B', alg], 
                        stdout = log_file, stderr = log_file, text = True)


def check_bam_type(bam):
    for read in pysam.AlignmentFile(bam, 'rb'):
        if read.is_paired:
            return 'BAMPE'
        else:
            return 'BAM'
        break


def check_alg_type(antibody):
    if re.findall(r'H[3-4]K', str(antibody)) or re.findall(
        r'H2[A-B]K', str(antibody)):
        return "--broad"
    else:
        return ""


def update_metadata(metadata, args):
    peak_files = []
    for sample in metadata.itertuples():
        peak_files.append(''.join([os.path.splitext(sample[5])[0], 
                                   '_peaks.xls']))
    metadata['Peaks'] = peak_files
    metadata['PeakCaller'] = 'macs'
    try:
        metadata.to_csv(args.metadata, sep='\t', index=False, )
    except:
        outfile = ''.join([os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '/examples/example_metadata.txt'])
        metadata.to_csv(outfile, sep = '\t')


def main():
    args = parse_args(sys.argv[1:])
    metadata, genome_dir, fig_dir, log_dir = input_params(args)
    chrom_sizes = genome_files(genome_dir)
    print('Performing peak calling')
    print('----------------------------------------------------------------')
    for sample in metadata.itertuples():
        peak_calling(sample, str(chrom_sizes[1].sum()), log_dir)
    if len(metadata.columns) == 7:
        update_metadata(metadata, args)


if __name__ == '__main__':
    main()
