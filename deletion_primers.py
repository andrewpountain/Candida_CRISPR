#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: andrewpountain

deletion_primers.py

Generates a file of knockout construct primers each with n bp overlap upstream and downstream of the gene
of interest.

Run:
python deletion_primers.py --sequences <sequence file in fasta format> --forward <sequence> --reverse <sequence>

Inputs:
--sequences: Sequence file containing coding sequences with 1 kb upstream and downstream flanking regions (can be .gz)
--forward: The forward binding site primer
--reverse: The reverse binding site primer
--overlap: The length of overlap with flanking regions to include
--flanks: The length of the sequences flanking the coding sequence in the input file

Outputs:
To standard out, the following tab-separated fields:
1. gene: Gene name
2. sequence_left: Sequence of the forward primer
3. sequence_right: Sequence of the reverse primer
4. issues: None if there are no issues, but if the flanking sequences contain non-ATCG bases it will specify which ones.

"""

import gzip
import argparse
import sys


def revcomp(x):
    """
    Reverse-complements DNA
    :param x: input sequence
    :return: reverse-complemented sequence
    """
    revdict = {"A": "T", "T": "A", "C": "G", "G": "C", "a": "t", "t": "a", "c": "g", "g": "c"}
    return "".join([revdict[i] for i in x[::-1]])

def checkATCG(x):
    """
    Test to see whether all letters in a sequence are A, T, C, or G.
    :param x: Input sequence
    :return: True if all letters in the input sequence are A, T, C, or G, otherwise False
    """
    legal_bases = {"A", "T", "C", "G", "a", "t", "c", "g"}
    return True if all(i in legal_bases for i in x) else False


# Parse the arguments
parser = argparse.ArgumentParser(description='Input details for indexing.')
parser.add_argument('--sequences', metavar='-s', type=str,
                    help='Provide the file name for the sequence fasta file (can use .gz files).')
parser.add_argument('--forward', metavar='-f', type=str,
                    help='Provide the sequence of the binding region of the forward primer.')
parser.add_argument('--reverse', metavar='-r', type=str,
                    help='Provide the sequence of the binding region of the reverse primer.')
parser.add_argument('--overlap', metavar='-n', type=int, default=50,
                    help='Provide the length of the overlap with the genomic DNA to include in the primer.')
parser.add_argument('--flanks', metavar='-l', type=int, default=1000,
                    help='Provide the length of the flanking regions upstream and downstream of the CDS.')

args = parser.parse_args()

# Open the sequences input file
file = gzip.open(args.sequences, "r") if args.sequences.endswith('.gz') else open(args.sequences)
print("gene\tsequence_left\tsequence_right\tissues")

# Now start to loop through, capturing sequences and printing primers as you go.
# Requires fasta format, so lines starting with ">" are used for the gene name, and subsequent lines are joined
# to form the sequence. When the next ">" header line is hit, it designs the primers and resets.
sequence = ''
gene = None

def output_primers(gene, sequence, args):
    upstream = sequence[args.flanks - args.overlap:args.flanks]
    downstream = sequence[len(sequence) - args.flanks:(len(sequence) - args.flanks) + args.overlap]

    if checkATCG(upstream) and checkATCG(downstream):
        # Now that you have the flanking regions from your gene, add on the forward and reverse binding
        # regions of the primer.
        forward_primer = upstream + args.forward
        reverse_primer = revcomp(downstream) + args.reverse
        print("\t".join([gene, forward_primer, reverse_primer, "None"]))
    else:
        issues = []  # Here any issues with non-ATCG bases will be added.
        if not checkATCG(upstream):
            issues.append("upstream sequence %s contains non-ATCG bases" % upstream)
        if not checkATCG(downstream):
            issues.append("downstream sequence %s contains non-ATCG bases" % downstream)
        print("\t".join([gene, "NA", "NA", ", ".join(issues)]))

for line in file.readlines():
    if line.startswith(b'>'):
        if len(sequence) > 0:
            output_primers(gene, sequence, args)
        sequence = ''
        gene = line.lstrip(b'>').split()[0].decode("utf-8")
    else:
        sequence += line.decode("utf-8").rstrip()
output_primers(gene, sequence, args)

file.close()
sys.exit()
