#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: andrewpountain

internal_primers.py

Designs internal primers for all coding sequences in the genome.
It makes use of the python interface with primer3-py, so that optimization is done by primer3.
The user supplies a config file with global input tags as described (http://primer3.org/manual.html).

Run: python3 internal_primers.py <sequences_fasta> <config_file> <gff_file> <number_primers>

Inputs:
sequences_fasta: fasta file (can be gzipped) containing sequences to be used as templates
config_file: file containing tab-separated key-value pairs with settings for global input tags
gff_file: file in gff format used to build exon structures of genes (the script is designed to look for only exonic binding sites)
number_primers: optional, defaults to 0, integer giving the number of primers to return

Outputs:
To stdout, a tab-separated table with the following fields (where multiple primers are returned, fields are comma-separated lists):
    1. gene: gene ID
    2. product_size: size of PCR amplicon(s)
    3. product_coordinates: hyphen-separated coordinates of PCR amplicon(s)
    4. sequence_left: DNA sequence of the forward primer(s)
    5. sequence_right: DNA sequence of the reverse primer(s)
    6. TM_left: Melting temperature of the forward primer(s)
    7. TM_right: Melting temperature of the reverse primer(s)
    8. GC_percent_left: GC percentage of the forward primer(s)
    9. GC_percent_right: GC percentage of the reverse primer(s)
    10. any_th_left: Tendency of forward primer(s) to bind to itself (see http://primer3.org/manual.html#PRIMER_MAX_SELF_ANY_TH)
    11. any_th_right: Tendency of reverse primer(s) to bind to itself (see http://primer3.org/manual.html#PRIMER_MAX_SELF_ANY_TH)
    12. end_th_left: Tendency of forward primer(s) to bind to itself at 3' end (see http://primer3.org/manual.html#PRIMER_MAX_SELF_END_TH)
    13. end_th_right: Tendency of reverse primer(s) to bind to itself at 3' end (see http://primer3.org/manual.html#PRIMER_MAX_SELF_END_TH)
    14. hairpin_th_left: Tendency of forward primer(s) to form a hairpin (see http://primer3.org/manual.html#PRIMER_MAX_HAIRPIN_TH)
    15. hairpin_th_right: Tendency of reverse primer(s) to form a hairpin (see http://primer3.org/manual.html#PRIMER_MAX_HAIRPIN_TH)
    16. other_allele: True/False for whether both alleles* can serve as a template for PCR
    17. fail: If primer design failes, this has a message saying this, otherwise "None"
    
*Note that in Assembly 22 of the Candida albicans SC5314 genome, there are two alleles for each coding sequence, 
denoted with suffixes "_A" and _B".    

"""
import primer3
import re
import sys
import gzip
from collections import OrderedDict


class Primer:
    """
    Class of objects representing individual primers
    """
    def __init__(self, sequence, TM, GC_percent, any_th, end_th, hairpin_th, direction):
        self.sequence = sequence
        self.TM = TM
        self.GC_percent = GC_percent
        self.any_th = any_th
        self.end_th = end_th
        self.hairpin_th = hairpin_th
        if direction in ["left","right"]:
            self.direction = direction
        else:
            raise ValueError("Primer direction must be either left or right.")

    def __len__(self):
        return len(self.sequence)

    def find_start(self, template):
        if getattr(self, "direction") == "left":
            position = template.find(self.sequence)
            return position + 1 if position != -1 else None
        elif getattr(self, "direction") == "right":
            position = template.find(revcomp(self.sequence))
            return position + len(self.sequence) if position != -1 else None

class PrimerPair:
    """
    Class representing a pair of objects of class Primer
    """
    def __init__(self, left, right, product_size):
        assert isinstance(left, Primer)
        assert isinstance(right, Primer)
        if getattr(left, "direction") == getattr(right, "direction"):
            raise ValueError("Primers in a primer pair cannot be in the same orientation.")
        self.left = left
        self.right = right
        self.product_size = int(product_size)

    def get_coords(self, template):
        """
        :return: Product coordinates i.e. the first coordinate of the left primer and the last coordinate of the
        right primer.
        """
        return self.left.find_start(template), self.right.find_start(template)

    def get_left(self):
        return self.left

    def get_right(self):
        return self.right

    def is_template(self, template):
        """
        Determines if a template matches or not.
        :param template: Input sequence
        :return: True if there is a match, False if not.
        """
        coords = self.get_coords(template)
        return True if all(coords) and coords[0] < coords[1] else False


class GeneStructure:
    """
    Class to define exon coordinates within a gene.
    Exons are defined as a tuple of coordinates relative to gene position.
    """
    def __init__(self, gene_id, start, end, sense):
        self.gene_id = gene_id
        self.start = int(start)
        self.end = int(end)
        self.sense = sense
        self.exons = []

    def get_exons(self):
        return self.exons

    def add_exon(self, exon_start, exon_end):
        """
        Produces an exon when supplied with the exon start and end position within the chromosome
        :param exon_start: start position within the chromosome
        :param exon_end: end position within the chromosome
        :return: Adds an exon to self.exons, a tuple containing the start position within the gene and the length.
        """
        exon_start = int(exon_start)
        exon_end = int(exon_end)

        exon_length = abs(exon_start - exon_end) # Primer3 seems to not count the first position in length, so the
                                                 # total length (+1) throws an error for being too long.

        if self.sense == "+":
            self.exons.append([exon_start - self.start + 1, exon_length])
        else:
            self.exons.append([self.end - exon_end + 1, exon_length])

    def get_region_list(self):
        """
        :return: list of regions formatted for SEQUENCE_PRIMER_PAIR_OK_REGION_LIST
        i.e. left_start, left_length, right_start, right_length
        """
        return [exon + exon for exon in self.exons]

    def __str__(self):
        return "%s: %s-%s %s, exons:%s" % (self.gene_id, self.start, self.end, self.sense, ",".join(str(exon) for exon in self.exons))


def revcomp(x):
    revdict = {"A": "T", "T": "A", "C": "G", "G": "C", "a": "t", "t": "a", "c": "g", "g": "c"}
    return "".join([revdict[i] for i in x[::-1]])


def import_config(file_name):
    return_dict = {}
    with open(file_name, "r") as file:
        for line in file.readlines():
            fields = line.rstrip().split("\t")
            if fields[0] == "PRIMER_PRODUCT_SIZE_RANGE":
                ranges = fields[1].split(' ')
                return_dict[fields[0]] = []
                for size_range in ranges:
                    return_dict[fields[0]].append([int(i) for i in size_range.split('-')])
            else:
                if re.match('\d+\.\d+', fields[1]):
                    return_dict[fields[0]] = float(fields[1])
                elif re.match('\d+', fields[1]):
                    return_dict[fields[0]] = int(fields[1])
                else:
                    return_dict[fields[0]] = fields[1]
    return return_dict


def find_primers(sequence_id, sequence_template, global_args, gene_coords_dict):
    """
    Runs the primer3.bindings.designPrimers function.
    :param sequence_id: Name of the sequence e.g. gene name
    :param sequence_template: The sequence of the template itself
    :param global_args: Dictionary of global arguments
    :return: A tuple of a boolean (whether or not successful) and a primer3 results object (dict).
    """
    seq_args = {'SEQUENCE_ID': sequence_id, 'SEQUENCE_TEMPLATE': sequence_template,
                'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': gene_coords_dict[sequence_id].get_region_list()}
    try:
        primer3_result = primer3.bindings.designPrimers(seq_args, global_args)
        return True, primer3_result
    except OSError as e:
        return False, "Primer design failed for gene %s. %s." % (sequence_id, str(e))


def parse_primers(primer3_result, number_primers):
    """
    Takes the primer results object derived from find_primers and returns a PrimerPair object.
    :param primer3_result: Output of find_primers function
    :param number_primers: Integer representing the number of primers you wish to return
    :return: Object of class PrimerPair
    """
    def make_primer(str, result, direction):
        """
        Helper function to initialise a Primer object
        :param str: Prefix to each argument
        :param result: primer3_result object
        :param direction: direction the primer is facing
        :return: Object of class Primer
        """
        try:
            return Primer(result[str + "_SEQUENCE"],
                          result[str + "_TM"],
                          result[str + "_GC_PERCENT"],
                          result[str + "_SELF_ANY_TH"],
                          result[str + "_SELF_END_TH"],
                          result[str + "_HAIRPIN_TH"],
                          direction)
        except KeyError:
            print(result)
            sys.exit()


    pair_list = []
    for i in range(number_primers):
        left_primer = make_primer("PRIMER_LEFT_" + str(i), primer3_result, "left")
        right_primer = make_primer("PRIMER_RIGHT_" + str(i), primer3_result, "right")
        pair_list.append(PrimerPair(left_primer, right_primer, primer3_result["PRIMER_PAIR_" + str(i) + "_PRODUCT_SIZE"]))
    return pair_list


def get_gene_structures(gff_file, info_tag = "ID", transcript_suffix = "-T"):

    gene_str_dict = {}
    file = gzip.open(gff_file, "r") if gff_file.endswith('.gz') else open(gff_file)
    for line in file.readlines():
        if not line.startswith("#"):
            fields = line.rstrip().split("\t")
            if fields[2] in ["gene", "pseudogene"]:
                info_fields = fields[8].split(";")
                for info in info_fields:
                    if info.startswith(info_tag + "="):
                        gene_id = info.split("=")[1]
                        break
                gene_str_dict[gene_id] = GeneStructure(gene_id, fields[3], fields[4], fields[6])

            elif fields[2] == "exon":
                info_fields = fields[8].split(";")
                for info in info_fields:
                    parent = None
                    if info.startswith("Parent="):
                        parent = info.split("=")[1]
                        break
                if parent:
                    gene_str_dict[parent.rstrip(transcript_suffix)].add_exon(fields[3], fields[4])
                else:
                    raise ValueError("The info tag 'Parent' is not found for exons in this gff file.")
    return gene_str_dict


def hits_both_alleles(primer_pair, gene, sequence_dict):
    """
    This function checks whether a primer pair hits the other allele from which it was derived.
    It works simply by switching the suffix from "_A" to "_B" and vice versa.
    If a gene name ends with neither, it assumes there is only one allele, so returns True.

    :param primer_pair: Object of class primer pair
    :param gene: Gene that the primer pair targets
    :return: True if it hits both, otherwise False.
    """
    if gene.endswith("_A"):
        other_allele = re.sub("_A$", "_B", gene)
    elif gene.endswith("_B"):
        other_allele = re.sub("_B$", "_A", gene)
    else:
        return True
    if other_allele in sequence_dict:
        return True if primer_pair.is_template(sequence_dict[other_allele]) else False
    else:
        return True


if __name__ == "__main__":
    sequences_file = sys.argv[1]
    config_file = sys.argv[2]
    gff_file = sys.argv[3]
    if len(sys.argv) > 4:
        number_primers = int(sys.argv[4])
    else:
        number_primers = 1

    # Import arguments from the config file
    global_args = import_config("/Users/andrewpountain/Documents/Lorenz_Postdoc/Knockout_library_generation/primer_design.config")

    # Build gene exon coordinates dictionary
    gene_coords_dict = get_gene_structures(gff_file)

    if global_args['PRIMER_NUM_RETURN']:
        if number_primers > global_args['PRIMER_NUM_RETURN']:
            raise ValueError("Number of primers cannot be greater than the 'PRIMER_NUM_RETURN' flag used in config.")

    first_line = "gene\tproduct_size\tproduct_coordinates\tsequence_left\tsequence_right\t" \
                 "TM_left\tTM_right\tGC_percent_left\tGC_percent_right\tany_th_left\tany_th_right\t" \
                 "end_th_left\tend_th_right\thairpin_th_left\thairpin_th_right\tother_allele\tfail"

    # Start looping through the fasta file, creating a sequence dict retaining the original order.
    file = gzip.open(sequences_file, "r") if sequences_file.endswith('.gz') else open(sequences_file)
    sequence = ''
    gene = None
    sequence_dict = OrderedDict()
    for line in file.readlines():
        if line.startswith(b'>'):
            if len(sequence) > 0:
                sequence_dict[gene] = sequence
            gene = line.lstrip(b'>').split()[0].decode("utf-8")
            sequence = ''
        else:
            sequence += line.decode("utf-8").rstrip()
    sequence_dict[gene] = sequence

    # Start looping through the sequence dict, performing primer3 as you go.
    print(first_line)
    for gene in sequence_dict.keys():
        sequence = sequence_dict[gene]
        # Run the primer find:
        successful, primer_result = find_primers(gene, sequence, global_args, gene_coords_dict)
        # Determine how many you want to retrieve (it can't be more than the actual number of primers returned.
        if successful:
            number_returned = primer_result['PRIMER_LEFT_NUM_RETURNED']
            if number_returned > 0:
                number_for_gene = number_primers if number_primers < number_returned else number_returned
            else:
                number_for_gene = 0

            primer_pairs = parse_primers(primer_result, number_for_gene)

            # Now I need to convert this into the line to be printed:
            print_line = [gene]
            # Product sizes
            print_line.append(",".join(str(getattr(pair, "product_size")) for pair in primer_pairs))
            # Product coordinates within the gene
            print_line.append(",".join("-".join(str(coord) for coord in pair.get_coords(sequence)) for pair in primer_pairs))
            for field in ["sequence", "TM", "GC_percent", "any_th", "end_th", "hairpin_th"]: # The primer-specific attributes you need
                print_line.append(",".join(str(getattr(pair.get_left(), field)) for pair in primer_pairs)) # Left primer attribute
                print_line.append(",".join(str(getattr(pair.get_right(), field)) for pair in primer_pairs)) # Right primer attribute
            print_line.append(",".join(str(hits_both_alleles(pair, gene, sequence_dict)) for pair in primer_pairs)) # Statement of whether it hits both alleles or not
            print_line.append("None")
            print("\t".join(print_line))
        else:
            print("\t".join([gene] + ['NA'] * 15 + [primer_result]))

    sys.exit()
