#!/usr/bin/env python3

"""
Created: August 23, 2023
Updated: August 23, 2023
Author(s): Todd Lenz, tlenz001@ucr.edu

Generates the ChIP-seq metadata file required as input for the DiffBind
R package. The output may be incorrect if replicates do not share a
similar naming convention. This script is meant to be used as a
component of the ChIPPy pipeline, but can be used independently if the
sorted BAM file names contain the string '_sorted.bam' and the XLS peak
files contain the string '_peaks.xls'.
"""


import sys
import argparse
import os
import pandas as pd


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',
                        '--input',
                        dest='input',
                        help='Directory containing sample folders.',
                        required=True)
    parser.add_argument('-t',
                        '--treatment',
                        dest='treatment',
                        help='Name or value indicating targeted antibody. '
                             'Targeted sample file names should contain this '
                             'string.',
                        required=True)
    parser.add_argument('-c',
                        '--control',
                        dest='control',
                        help='Name or value indicating control antibody or '
                             'input samples. Control file names should '
                             'contain this string.',
                        required=True)


def get_file_paths(input, bam_string):
    """ Returns paths to files in 'input' directory. """
    file_paths = []
    for root, _, files in os.walk(input):
        for file in files:
            if file.endswith(bam_string):
                full_path = os.path.join(root, file)
                file_paths.append(full_path)
    return file_paths


def get_substrings(samples, factor):
    """ 
    Removes underscores and dashes from the ends and any duplicated 
    underscores and dashes with single copies from sample names.
    """
    samples = [s.replace(factor, '').strip('_').strip('-') for s in samples]
    while any('__' in s for s in samples):
        samples = [s.replace('__', '_') for s in samples]
    while any('--' in s for s in samples):
        samples = [s.replace('--', '-') for s in samples]
    return samples


def split_samples(samples):
    """ Groups sample conditions based on sample names. """
    for sample in samples:
        for l in range(len(sample), 0, -1):
            for s in range(len(sample) - l + 1):
                substr = sample[s:s + l]
                sub = [cs for cs in samples if substr in cs]
                if len(sub) > 1:
                    samples = [x for x in samples if x not in sub]
                    return samples, substr, sub


def main():
    args = parse_args(sys.argv[1:])

    file_list = [os.path.abspath(f) for f in os.listdir(args.input) if 
                 os.path.isdir(os.path.join(args.input, f))]
    sampleIDs = [os.path.basename(f) for f in file_list if args.treatment in f]
    controlIDs = [os.path.basename(f) for f in file_list if args.control in f]
    bam_files = get_file_paths(args.input, '_sorted.bam')
    bamReads = [f for f in bam_files if args.treatment in f]
    bamControls = [f for f in bam_files if args.control in f]
    Peaks = get_file_paths(args.input, '_peaks.xls')

    samples = get_substrings(sampleIDs, args.treatment)
    conditions = {}
    while len(samples) > 0:
        samples, cond, subsamples= split_samples(samples)
        conditions[cond] = subsamples
    cdf = pd.DataFrame(conditions)
    cdf = pd.melt(cdf)

    df = pd.DataFrame({'SampleID': sampleIDs, 
                       'Factor': args.treatment, 
                       'Condition': cdf['variable'], 
                       'Replicate': cdf.groupby('variable').cumcount() + 1, 
                       'bamReads': bamReads, 
                       'ControlID': controlIDs, 
                       'bamControl': bamControls, 
                       'Peaks': Peaks, 
                       'PeakCaller': 'macs'})
    df.to_csv("chip_metadata.txt", sep='\t', )


if __name__ == '__main__':
    main()
