#!/usr/bin/env python
"""
Created June 13, 2022
Updated August 8, 2023

@author: Todd Lenz, tlenz001@ucr.edu

Plots genomewide coverage for ChIP-seq data. Track lengths are 
normalized with respect to the longest chromosome--i.e. the longest 
chromosome will fill the width of the figure and all other chromosomes 
are proportionally plotted against it. Control reads (input/IGG) should 
already be subtracted from the sample data.

Required input arguments:
1. bed--BED file(s) containing genome-wide per-base coverage. 
   To generage this file from a BAM file:
   bedtools genomecov -d -ibam ${i}_sorted.bam > ${i}_sorted.bed

Optional input arguments:
1. gff--Standard GFF file for genome/organism of interest.
   (https://www.ensembl.org/info/website/upload/gff.html?redirect=no)
2. samples--List of sample names to annotate the y-axes of barplots. 
   List should be the sample length as the list of BED files and plots 
   will be annotated in the same order as BED files. If no sample 
   names are given, the BED file names are used without the extension.
3. output--Directory and name for output .pdf file. If no output is 
   given, the output file is named 'output.pdf' and saved in the 
   directory of the first BED input.
4. resolution--Integer indicating the binning resolution (default = 10).
5. centromeres--Tab-delimited file containing centromere coordinates 
   for each chromosome. File format:
   chromosome_name    start_coordinate    end_coordinate
   Chromosome names ('chromosome_name') should match the names in the 
   BED files.
6. gene_list--List of genes to plot, separated by spaces, commas, 
   newlines, or tabs. Genes will be extracted from the provided GFF 
   file.
7. normalization--To normalize data in BED files so that samples can be 
   directly compared, indicate a method for normalization. Data can be 
   normalized using counts-per-million (CPM). More methods will be 
   added in the future.
8. ymax--If argument is provided, all chromosome plots will use the 
   maximum y-value for the entire binned dataset. If not provided, each 
   chromosome will use the maximum value for that chromosome.
"""


import sys
import os
import math
import re
import argparse
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-b',
                        '--bed',
                        dest='bed',
                        help='Tab-delimited BED files with per-base '
                             'read coverage (bedtools genomecov -d).',
                        nargs='+',
                        required=True)
    parser.add_argument('-g',
                        '--gff',
                        dest='gff',
                        help='Gene data in GFF format. Coding regions are '
                             'used to determine which coordinates within BED '
                             'files to plot.',
                        required=False,
                        default=None)
    parser.add_argument('-s',
                        '--samples',
                        dest='samples',
                        help='Sample name(s) used to annotate barplot(s).',
                        nargs='+',
                        required=False,
                        default=None)
    parser.add_argument('-o',
                        '--output',
                        dest='output',
                        help='Output file path and name. If no file is given, '
                             'the output will be saved in the same directory '
                             'as the first BED file with a default file name.',
                        required=False,
                        default=None)
    parser.add_argument('-r',
                        '--resolution',
                        dest='resolution',
                        help='Resolution at which to bin read count data. '
                             'Binning will begin at coordinate 1 and continue '
                             'to the end of the gene, therefore the last bin '
                             'may be less than the given value.',
                        required=False,
                        default=10000,
                        type=int)
    parser.add_argument('-c',
                        '--centromeres',
                        dest='centromeres',
                        help='List of centromere coordinates.',
                        required=False,
                        default=None)
    parser.add_argument('-l',
                        '--gene_list',
                        dest='gene_list',
                        help='List of genes to plot, separated by spaces, '
                             'commas, newlines or tabs. Genes will be '
                             'extracted from the provided GFF file.',
                        required=False,
                        default=None)
    parser.add_argument('-n',
                        '--normalize',
                        dest='normalize',
                        help='Normalize read counts for each sample to '
                             'account for differences in sequencing depth. '
                             'Data can be normalized using CPM (counts-per-'
                             'million) normalization.',
                        required=False,
                        choices=['CPM'],
                        default=None)
    parser.add_argument('-y',
                        '--ymax',
                        dest='ymax',
                        help='Set equal maximum y-value for all chromosomes.',
                        required=False,
                        action='store_false',
                        default=None)
    return parser.parse_args()


def input_params(args):
    """
    Parse input arguments, loading bed files into dataframes. If 
    provided, the samples, resolution, and distance arguments will also 
    be parsed and genes will be extracted from an input gff file.
    """
    genes, centromeres = None, None
    if args.gff:
        genes = filter_gff(pd.read_csv(args.gff, sep='\t', header=None, 
                                       usecols=[0,2,3,4,8]))
        if args.gene_list:
            genes = extract_genes(genes, args.gene_list)
        genes.columns = range(len(genes.columns))
    if args.centromeres:
        centromeres = pd.read_csv(args.centromeres, sep='\t', header=None)
    samples = []
    for i in args.bed:
        bed = pd.read_csv(i, sep='\t', header=None)
        try:
            df = df.merge(bed, on=[0,1])
        except:
            df = bed
            if not args.output:
                out = ''.join([os.path.dirname(i), '/ChIP_barplot.pdf'])
            else:
                out = args.output
        if not args.samples or len(args.samples) != len(args.bed):
            samples.append(os.path.basename(i)[:-4])
    df.columns = range(len(df.columns))
    if len(samples) == 0:
        samples = args.samples
    if args.normalize:
        df = normalize_df(df, args.normalize)
    res = int(args.resolution)
    return genes, centromeres, df, out, samples, res


def filter_gff(gff):
    """ Extracts gene accessions, names and descriptions. """
    genes = gff[(gff[2] == 'protein_coding_gene') | (gff[2] == 'ncRNA_gene') | 
                (gff[2] == 'pseudogene')].reset_index(drop=True).copy()
    genes[9] = genes[8].str.extract(r'ID=(.*?);', expand=True)
    genes[10] = genes[8].str.extract(r'Name=(.*?);', expand=True)
    genes[11] = genes[8].str.extract(r'description=(.*?);', expand=True)
    genes.drop([2,8], axis=1, inplace=True)
    genes[11] = genes[11].str.replace(r'%2C', ',', regex=True)
    return genes


def extract_genes(genes, gene_list):
    """
    Filters genes by list provided in text file. Gene names can be 
    delimited by spaces, commas, newlines or tabs.
    """
    delimiter = check_delimiter(gene_list)
    list_of_genes = []
    with open(gene_list, 'r') as f:
        for line in f:
            items = line.strip().split(delimiter)
            list_of_genes.extend(items)
    genes = genes[genes[9].isin(list_of_genes)]
    return genes


def check_delimiter(gene_list):
    """ Identify delimiter for list of genes. """
    with open(gene_list, 'r') as f:
        lines = [f.readline().strip() for _ in range(5)]
    space_count = sum(line.count(' ') for line in lines)
    comma_count = sum(line.count(',') for line in lines)
    newline_count = sum(line.count('\n') for line in lines)
    tab_count = sum(line.count('\t') for line in lines)
    max_count = max(space_count, comma_count, newline_count, tab_count)
    if space_count == max_count:
        return ' '
    if comma_count == max_count:
        return ','
    elif newline_count == max_count:
        return '\n'
    elif tab_count == max_count:
        return '\t'
    else:
        print('Can\'t determine delimiter of gene_list.')


def normalize_df(df, norm):
    """ Normalize data from BED files. """
    if norm == 'CPM':
        for i in range(2,len(df.columns)):
            mmr = df[i].sum() / 1000000
            df[i] = df[i].div(mmr)
    else:
        print('Provided method of normalization is not valid so data will not '
              'be normalized. Run script with -h flag to see valid '
              'normalization methods.')
        pass
    return df


def data_binning(df, res):
    """ Bins read counts for each chromosome. """
    df[len(df.columns)] = df.iloc[:, 1] // res + 1
    df = df[list(df.columns)[0:1] + 
            list(df.columns)[-1:] + 
            list(df.columns)[2:-1]]
    df.columns = range(len(df.columns))
    df = df.groupby([0,1], as_index=False)[list(df.columns)[2:]].sum()
    return df


def custom_ax_params(ax, sample, max_yval, max_xval):
    """ Set font parameters and axes labels, limits and tick marks. """
    font = {'family':'serif',
            'color':'black',
            'weight':'bold',
            'size':10}
    ax.set_ylim(bottom=0, top=max_yval)
    plt.ylabel(sample, rotation='horizontal', labelpad=20, 
               fontdict=font, ha='right')
    ax.set_xlim(left=0, right=max_xval)
    plt.tick_params(axis='x', bottom=False, labelbottom=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1)
    ax.spines['left'].set_capstyle('round')
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['bottom'].set_capstyle('round')
    ax.spines['bottom'].set_position('zero')
    return ax


def annotate_genes(ax, genes, res, max_yval):
    """ Annotate each plot with bars representing genes. """
    for g in genes.itertuples():
        ax.arrow(g[2] // res + 1, -max_yval * .2, (g[3]-g[2]) // res + 1, 
                 0, width=max_yval/4, head_width=0, head_length=0, 
                 facecolor='#FF0000', edgecolor='#FF0000', 
                 length_includes_head=True, clip_on=False)
    return ax


def main():
    args = parse_args(sys.argv[1:])
    genes, centromeres, df, out, samples, res = input_params(args)
    df = data_binning(df, res)
    if args.ymax is not None:
        max_yval = df[list(df.columns)[2:]].max().max()
    max_xval = max(df[0].value_counts())
    pdf = PdfPages(out)
    num_plots = len(df[0].unique())*(len(samples)+1)
    plot_range = [*range(num_plots)]
    plot_idx = 0
    sample_colors = ['#D81B60', '#1E88E5', '#FFC107']
    fig = plt.figure()
    fig.set_figheight(num_plots)
    fig.set_figwidth(20)
    for chr, idx in zip(df[0].unique(), range(len(df[0].unique()))):
        chr_df = df[df[0] == chr].copy()
        if args.ymax is None:
            max_yval = chr_df[list(chr_df.columns)[2:]].max().max()
        for i in [*range(len(samples))]:
            ax = plt.subplot2grid((num_plots, 20), (plot_range[plot_idx], 0), 
                                  colspan=int(math.ceil(max(chr_df[1]) / 
                                                        max_xval * 20)), 
                                                        rowspan=1)
            barplt = plt.bar(np.array(chr_df[1]), np.array(chr_df[2+i]), 
                             width=1, color=sample_colors[i])
            ax = custom_ax_params(ax, samples[i], max_yval, max(chr_df[1]))
            plot_idx += 1
            if i == len(samples) - 1:
                plot_idx += 1
                if genes is not None:
                    ax = annotate_genes(ax, genes[genes[0] == chr_df[0].iloc[0]
                                                  ].reset_index(drop=True), 
                                                  res, max_yval)
    fig.subplots_adjust(hspace=0)
    fig.tight_layout(h_pad=0)
    pdf.savefig()
    plt.close()
    pdf.close()


if __name__ == '__main__':
    main()
