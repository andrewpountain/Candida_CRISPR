#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: andrewpountain

make_primer_files.py

Script that takes the outputs of three separate scripts and turns these into formatted text files for each gene.

python make_primer_files.py <sgRNA_primers.py output> <deletion_primers.py output> <internal_primers.py output> <output directory>
"""

import sys
import os
import pandas as pd

sgRNA_df = pd.read_csv(sys.argv[1], sep="\t")
deletion_df = pd.read_csv(sys.argv[2], sep="\t")
internal_df = pd.read_csv(sys.argv[3], sep="\t")
output_dir = os.path.abspath(sys.argv[4])

genes = list(sgRNA_df["gene"])



for gene in genes:
    # Specify an output file
    outfile = open(output_dir + "/" + gene + "_CRISPR_primers.txt", "w")


    # Guide RNA
    sgRNA_target = sgRNA_df.loc[sgRNA_df["gene"] == gene, "target"].iloc[0]
    pSNR52_R = sgRNA_df.loc[sgRNA_df["gene"] == gene, "pSNR52_R"].iloc[0]
    scaffold_F = sgRNA_df.loc[sgRNA_df["gene"] == gene, "scaffold_F"].iloc[0]
    target_sense = sgRNA_df.loc[sgRNA_df["gene"] == gene, "target_sense"].iloc[0]
    sense_vs_gene = sgRNA_df.loc[sgRNA_df["gene"] == gene, "sense_vs_gene"].iloc[0]

    outfile.write("Gene: %s\n"
          "\n"
          "Guide RNA:\n"
          "sgRNA target sequence: %s, (%s) sense compared to genome, (%s) sense compared to gene\n"
          "pSNR52_R primer: %s\n"
          "scaffold_F primer: %s\n\n" % (gene, sgRNA_target, target_sense, sense_vs_gene, pSNR52_R, scaffold_F))


    # Deletion primers
    outfile.write("Deletion construct primers:\n")
    if gene in deletion_df["gene"].values:
        deletion_issues = deletion_df.loc[deletion_df["gene"] == gene, "issues"].iloc[0]
        if deletion_issues != "None":
            outfile.write(deletion_issues + "\n")
        else:
            upstream = deletion_df.loc[deletion_df["gene"] == gene, "sequence_left"].iloc[0]
            downstream = deletion_df.loc[deletion_df["gene"] == gene, "sequence_right"].iloc[0]
            outfile.write("Forward primer: %s\nReverse primer: %s\n\n" % (upstream, downstream))
    else:
        outfile.write("No deletion primers found for this gene.\n\n")

    # Internal primers
    if gene in internal_df["gene"].values:
        fail_description = internal_df.loc[internal_df["gene"] == gene, "fail"].iloc[0]
        if fail_description == "None":
            product_size_list = internal_df.loc[internal_df["gene"] == gene, "product_size"].iloc[0].split(",")
            product_coordinates_list = internal_df.loc[internal_df["gene"] == gene, "product_coordinates"].iloc[0].split(",")
            left_list = internal_df.loc[internal_df["gene"] == gene, "sequence_left"].iloc[0].split(",")
            right_list = internal_df.loc[internal_df["gene"] == gene, "sequence_right"].iloc[0].split(",")
            TM_left_list = internal_df.loc[internal_df["gene"] == gene, "TM_left"].iloc[0].split(",")
            TM_right_list = internal_df.loc[internal_df["gene"] == gene, "TM_right"].iloc[0].split(",")
            GC_percent_left_list = internal_df.loc[internal_df["gene"] == gene, "GC_percent_left"].iloc[0].split(",")
            GC_percent_right_list = internal_df.loc[internal_df["gene"] == gene, "GC_percent_right"].iloc[0].split(",")
            other_allele_list = internal_df.loc[internal_df["gene"] == gene, "other_allele"].iloc[0].split(",")

            primers_found = False
            for i in range(len(product_size_list)):
                if bool(other_allele_list[i]):
                    if not primers_found:
                        outfile.write("Internal primers:\n\nPair number\tProduct size\tProduct coordinates\tForward sequence\tReverse sequence\tForward TM\tReverse TM\tForward GC\tReverse GC\n")
                    primers_found = True
                    outfile.write("\t".join(["Pair " + str(i + 1), product_size_list[i], product_coordinates_list[i],
                        left_list[i], right_list[i], str(round(float(TM_left_list[i]), 2)),
                        str(round(float(TM_right_list[i]), 2)), str(round(float(GC_percent_left_list[i]), 1)),
                        str(round(float(GC_percent_right_list[i]), 1))]) + "\n")
            if not primers_found:
                outfile.write("No internal primers found for this gene. Check the original table, internal_primers.txt, to determine whether this is an issue of allele specificity.\n")
        else:
            outfile.write("No internal primers found for this gene. " + fail_description + "\n")
    else:
        outfile.write("No internal primers found for this gene. It may have been missing from the original fasta file used to design primer sequences.")
    outfile.close()

sys.exit()