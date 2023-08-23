#!/usr/bin/env python3

import re, os, sys, csv
import pandas as pd
import numpy as np


# Run script with two arguments:
# 1. GFF file containing gene info for organism of interest
# 2. MACS2/3 narrowPeak or broadPeak file containing peak calling data


def assign_args(args):

	# Iterate over input arguments
	for i in args:

		# Input is gff containing gene info
		if os.path.splitext(i)[1] == ".gff":

			# Generate pandas dataframe from input gff
			gff = pd.read_csv(i, sep='\t', header=None)

			# Filter gff dataframe to only include protein coding genes
			genes = gff.loc[(gff[2] == "protein_coding_gene") | (gff[2] == "ncRNA_gene")]

			# Extract gene accession numbers
			accession = genes[8].str.extract(r'ID=(.*?);', expand=True)

			# Extract gene names
			name = genes[8].str.extract(r'Name=(.*?);', expand=True)

			# Extract gene descriptions
			description = genes[8].str.extract(r'description=(.*?);', expand=True)

			# Concatenate new columns to genes df
			genes = pd.concat([genes, accession, name, description], axis=1)

			# Re-index rows and column names
			genes.columns = range(genes.columns.size)
			genes.reset_index(drop=True, inplace=True)

			# Replace text in gene description
			genes[11] = genes[11].str.replace(r'%2C', ',', regex=True)
			genes[11] = genes[11].str.replace(r'+', ' ', regex=True)

			# Sort genes by chromosome and start position
			genes.sort_values([0,3], ascending=True, inplace=True)
			genes.reset_index(drop=True, inplace=True)

		# Input is broadPeak or narrowPeak file containing peak info
		elif os.path.splitext(i)[1] == ".narrowPeak" or os.path.splitext(i)[1] == ".broadPeak":

			# Generate pandas dataframe from input peaks file
			peaks = pd.read_csv(i, sep='\t', header=None)

			# Sort peaks by chromosome and start position
			peaks.sort_values([0, 1], ascending=True, inplace=True)
			peaks.reset_index(drop=True, inplace=True)

			# Assign output directory
			outdir = os.path.dirname(os.path.abspath(i))

			# Assign output file name
			outfile = "".join([os.path.splitext(os.path.basename(i))[0], "_mapped.txt"])

		elif os.path.splitext(i)[1] == ".bed" or os.path.splitext(i)[1] == ".txt":

			# Generate pandas dataframe from input peaks file
			peaks = pd.read_csv(i, sep='\t')
			peaks.columns = range(len(peaks.columns))

			# Sort peaks by chromosome and start position
			peaks.sort_values([0, 1], ascending=True, inplace=True)
			peaks.reset_index(drop=True, inplace=True)

			# Assign output directory
			outdir = os.path.dirname(os.path.abspath(i))

			# Assign output file name
			outfile = "".join([os.path.splitext(os.path.basename(i))[0], "_mapped.txt"])

		# If any input is not a GFF or narrowPeak/broadPeak format, 
		# an error is raised and the script exits immediately.
		else:
			print()
			print("".join(["Input file ", i, " is not a recognized input format."]))
			print("Run script with standard GFF and narrowPeak/broadPeak files.")
			print()
			sys.exit(0)

	return genes, peaks, outdir, outfile


def check_overlaps(sp, sg):
	overlap = 0

	# Create list of values covered by peak
	sp_vals = list(range(sp[2],sp[3]))

	# Create empty list to store info for gene with largest overlap of peak
	peak_overlap = []

	# Iterate over genes that overlap peak
	for g in sg.itertuples():

		# Create list of values covered by gene
		g_vals = list(range(g[4],g[5]))

		# Find values covered by both peak and gene, i.e. length of overlap
		ol = len(list(set(sp_vals) & set(g_vals)))

		# If current gene has largest overlap, store info in peak_overlap list
		if ol > overlap:
			overlap = ol
			peak_overlap = g[10:]
		else:
			continue
	return peak_overlap


def main():

	# Assign input and output args
	genes, peaks, outdir, outfile = assign_args(sys.argv[1:])

	# Generate empty list to store genes mapped to peaks
	gene_peaks = []

	# Iterate over chromosomes separately to speed up runtime
	for chrom in peaks[0].unique():

		# Subset peaks and genes by current chromosome
		sub_peaks = peaks[peaks[0] == chrom].copy()
		sub_genes = genes[genes[0] == chrom].copy()

		# Iterate over peaks in subset
		for sp in sub_peaks.itertuples():

			# Subset genes that overlap the current peak
			sg = sub_genes[(sub_genes[3] - 1000 < (sp[3] - sp[2]) / 2 + sp[2]) & (sub_genes[4] > (sp[3] - sp[2]) / 2 + sp[2])].copy()

			# If no genes overlap the current peak append list with "intergenic"
			if len(sg.index) == 0:
				gene_peaks.append(["intergenic", "", ""])

			# If only 1 gene overlaps the current peak append list with gene info
			elif len(sg.index) == 1:
				gene_peaks.append(sg.iloc[:, 9:].values.flatten().tolist())

			# If multiple genes overlap the current peak,
			# calculate gene with greatest overlap and append list with gene info
			else:
				gene_peaks.append(sg.iloc[:, 9:].values.flatten().tolist())
				#gene_peaks.append(check_overlaps(sp, sg))

	# Create pandas dataframe from list of genes overlapping peaks
	gene_peaks_df = pd.DataFrame(gene_peaks)

	# Concatenate gene df with peak df
	peaks_mapped = pd.concat([peaks, gene_peaks_df], axis=1)

	# Write genes mapped to peaks to tab-delimited output
	peaks_mapped.to_csv("".join([outdir,"/",outfile]), sep="\t", index=False, header=False)


if __name__ == "__main__":
	main()