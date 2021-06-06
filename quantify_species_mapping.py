#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: andrewpountain
quantify_species_mapping.py
Quantifies mapping numbers to each of a given number of sets of genome features.
This is when, for example, you have mixed species data (e.g. dual host-pathogen RNA-seq),
and you have aligned all the reads to a concatenated, multi-species genome. This will
give you stats on the number of reads aligned to each species in there.
Argument 1 is a comma-separated list of files, each of which should have a newline-
separated list of chromosomes (N.B. if a chromosome appears in multiple files, it will
be counted multiple times).
Subsequent arguments are listed bam files you want to quantify.
Bam files are loaded in and these are filtered as follows:
1. total: all alignment records are counted
2. paired: all read 1 alignment records that are paired, and primary or unmapped
3. unmapped: of paired, all that are unmapped
4. properly_paired: of paired, all that are properly paired
5. species_mapping: dictionary of the chromosome sets you provided in the first argument,
					with values equalling the number of records in properly_paired that
					map to features in that chromosome set
Note that technically only read 1 records are quantified for all fields except total, 
but as only paired reads are counted this equates to number of read pairs.
"""

import pysam
import sys

def read_chromosomes(file):
	# Turns supplied file into a list of features
	chromosomes_list = []
	with open(file, "r") as file:
		for line in file.readlines():
			if not line.isspace():
				chromosomes_list.append(line.rstrip())
	return chromosomes_list

chromosome_files = sys.argv[1].split(",")
chromosome_dict = {}
for file in chromosome_files:
	chromosome_dict[file] = read_chromosomes(file)
bam_files = sys.argv[2:]

header = ["file", "total", "paired", "unmapped", "properly_paired"] + chromosome_files
print("\t".join(header))

for bam_file in bam_files:
	bam = pysam.AlignmentFile(bam_file, "rb") # Note that right now this only works for bam-compressed alignment files, but easily modified for e.g. sam and cram
	total = 0
	paired = 0
	properly_paired = 0
	unmapped = 0
	species_mapping = {chromosomes : 0 for chromosomes in chromosome_dict}
	for read in bam:
		total += 1
		if read.is_read1 and read.is_paired and not read.is_secondary and not read.is_supplementary:
			paired += 1
			if read.is_unmapped:
				unmapped += 1
			elif read.is_proper_pair:
				properly_paired += 1
				for key in chromosome_dict.keys():
					if read.reference_name in chromosome_dict[key]:
						species_mapping[key] += 1
						
	print("%s\t%d\t%d\t%d\t%d\t%s" % (bam_file, total, paired, unmapped, properly_paired, "\t".join([str(species_mapping[chromosomes]) for chromosomes in chromosome_dict.keys()])))
		
