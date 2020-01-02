#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: andrewpountain

sgRNA_primers.py

Designs primers to amplify CRISPR guide RNA constructs.

The sgRNA is chosen based on the following prioritizing rules:
    1. The most central guide RNA within the gene is chosen
    2. If there is a tie, guides in the same orientation as the gene are preferred (for convenience)
    3. If the sense is the same, the choice is made arbitarily based on the order of testing.
    4. If a sequence contains non-ATCG bases, it is ignored.

Run: python sgRNA_primers.py <sgRNAs>

The <sgRNAs> file is a text file containing annotation for guide RNAs.

This script defines a class, GenomicFeature, that has the coordinates and sense of a feature. The key methods of this
class are "get_midpoint", which outputs the middle position of that gene, and "midpoint distance", which calculates
the number of bases between the midpoint of that feature and another. This is used to choose the sgRNA target that
is most central within a gene.

Two objects inherit:
    GeneCoords: The coordinates for the whole gene. Objects of class GuideRNA can be added to the targets attribute.
    GuideRNA: An additional method is cut_PAM (the NGG PAM sequence should not be included in the final primer).

Prints the following tab-separated fields:
1. gene: The name of the gene
2. target: The sequence of the sgRNA target sequence
3. midpoint_distance: The distance of the middle of the targets sequence from the middle of the gene
4. pSNR52_R: Sequence of the pSNR52 reverse primer
5. scaffold_F: Sequence of the scaffold forward primer
6. target_sense: Sense of the target sequence relative to the genome
7. sense_vs_gene: Sense of the target sequence relative to the gene

"""

import sys
import re


class GenomicFeature:
    """
    Objects of class feature contain key coordinate information, including chromosomal coordinates and sense.
    The method get_midpoint also reports the midpoint coordinate for the feature.
    """
    def __init__(self, chrom, start, end, sense):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        if sense in ("+", "-"):
            self.sense = sense
        else:
            raise ValueError("'%s' can not be used as a value for gene sense." % sense)

    def get_chrom(self):
        return self.chrom.rstrip("A|B")

    def get_coords(self):
        return (self.start, self.end)

    def get_sense(self):
        return self.sense

    def get_midpoint(self):
        return self.start + (self.end - self.start)//2

    def midpoint_distance(self, other):
        assert isinstance(other, GenomicFeature)
        if self.get_chrom() != other.get_chrom():
            raise ValueError("You can not compare the distances between midpoints for"
                             " features on separate chromosomes:\n%s\n%s" % (self.__str__(), other.__str__()))
        return abs(self.get_midpoint() - other.get_midpoint())


class GeneCoords(GenomicFeature):
    """
    Objects of class gene_coords contain key coordinate information including gene ID, chromosomal coordinates
    and sense. The method get_midpoint also reports the midpoint coordinate for the gene.
    """
    def __init__(self, gene_id, chrom, start, end, sense):
        GenomicFeature.__init__(self, chrom, start, end, sense)
        self.gene_id = gene_id
        self.targets = set()

    def get_name(self):
        return self.gene_id

    def add_target(self, target):
        assert isinstance(target, GuideRNA)
        self.targets.add(target)

    def get_targets(self):
        return self.targets

    def print_targets(self):
        for target in self.targets:
            print(target)

    def __str__(self):
        return "%s\t%s:%d-%d\t%s" % (self.gene_id, self.chrom, self.start, self.end, self.sense)


class GuideRNA(GenomicFeature):
    def __init__(self, sequence, chrom, start, end, sense):
        GenomicFeature.__init__(self, chrom, start, end, sense)
        self.sequence = sequence

    def get_sequence(self):
        return self.sequence

    def cut_PAM(self):
        """
        :return: Sequence with the NGG PAM sequence removed.
        """
        if re.search("[ATCG]GG$", self.sequence):
            return self.sequence[:len(self.sequence) - 3]
        else:
            raise ValueError("The sgRNA %s does not end in PAM sequence NGG." % self.sequence)

    def __str__(self):
        return "%s\t%s:%d-%d\t%s" % (self.sequence, self.chrom, self.start, self.end, self.sense)


def process_line(gene_dict, line, allele_specific = False):
    """
    This processes input lines into GuideRNA and GeneCoords objects.
    :param gene_dict: The main gene_dict dictionary from the script.
    :param line: The input line of text from the sgRNA file.
    :param allele_specific: If it is true, it will create separate GuideRNA and GeneCoords objects for each allele.
                            If false, only the _A allele is included. Note I haven't implemented the True case in this
                            script, as for now I am just trying to select guides that hit both.
    :return: None
    """
    fields = line.rstrip().split("\t")
    sequence = fields[1]
    sgRNA_chrom = tuple(fields[2].split(","))
    sgRNA_start = tuple(fields[3].split(":"))
    sgRNA_end = tuple(fields[4].split(":"))
    sgRNA_sense = tuple(fields[5].split(","))
    gene_chrom = tuple(fields[6].split(","))
    gene_start = tuple(fields[7].split(":"))
    gene_end = tuple(fields[8].split(":"))
    gene_sense = tuple(fields[9].split(","))
    gene_id = tuple(fields[11].split(","))
    
    if gene_id[0] not in gene_dict:
        gene_dict[gene_id[0]] = GeneCoords(gene_id[0], gene_chrom[0], gene_start[0], gene_end[0], gene_sense[0])
    gene_dict[gene_id[0]].add_target(GuideRNA(sequence, sgRNA_chrom[0], sgRNA_start[0], sgRNA_end[0], sgRNA_sense[0]))

    # Now to the second allele if specified:
    if allele_specific:
        if gene_id[1] not in gene_dict:
            gene_dict[gene_id[1]] = GeneCoords(gene_id[1], gene_chrom[1], gene_start[1], gene_end[1], gene_sense[1])
        gene_dict[gene_id[1]].add_target(GuideRNA(sequence, sgRNA_chrom[1], sgRNA_start[1], sgRNA_end[1], sgRNA_sense[1]))


def revcomp(x):
    revdict = {"A": "T", "T": "A", "C": "G", "G": "C", "a": "t", "t": "a", "c": "g", "g": "c"}
    return "".join([revdict[i] for i in x[::-1]])


def checkATCG(x):
    legal_bases = {"A", "T", "C", "G", "a", "t", "c", "g"}
    return True if all(i in legal_bases for i in x) else False


if __name__ == "__main__":
    # Define the generic parts of the primers:
    generic_pSNR52_R = "CAAATTAAAAATAGTTTACGCAAGTC"
    generic_scaffold_F = "GTTTTAGAGCTAGAAATAGCAAGTTAAA"

    # Write the header line for the output file:
    print("gene\ttarget\tmidpoint_distance\tpSNR52_R\tscaffold_F\ttarget_sense\tsense_vs_gene")

    # Start looping through lines and building up a gene dictionary, each containing a set of guide RNAs associated.
    gene_dict = {}
    with open(sys.argv[1]) as file:
        file.readline()
        for line in file.readlines():
            process_line(gene_dict, line)

    # Now run through each one and find the best PAM sequence.
    for gene in gene_dict.keys():
        gene_obj = gene_dict[gene]
        best_target = None

        for target in gene_obj.get_targets():
            if not checkATCG(target.sequence):
                continue  # This excludes reads that contain non-ATCG bases.
            try:
                distance = target.midpoint_distance(gene_obj)
            except ValueError as e: # I'll leave this in for now, but it was only an issue when writing this.x
                print(e)
                print(gene_obj)
                gene_obj.print_targets()
                sys.exit()
            if best_target:  # (i.e. if it's not the first target you've looked at)
                # This bit does the search for the best one
                if distance < best_distance:  # First on the basis of distance
                    best_target = target
                    best_distance = distance
                elif distance == best_distance:  # Then on the basis of sense
                    gene_sense = gene_obj.get_sense()
                    if target.get_sense() == gene_sense and best_target.get_sense() != gene_sense:
                        best_target = target
            else:  # The case for the first target looked at (by definition the best looked at so far).
                best_target = target
                best_distance = distance
        if best_target:  # (i.e. if you've found an sgRNA matching criteria)
            # Add the sgRNA sequence (in the right orientation) to the generic sequences, removing the PAM sequence.
            pSNR52_R = revcomp(best_target.cut_PAM()) + generic_pSNR52_R
            scaffold_F = best_target.cut_PAM() + generic_scaffold_F
            if gene_obj.get_sense() == best_target.get_sense():
                sense_vs_gene = "+"
            else:
                sense_vs_gene = "-"
            print("\t".join([gene, best_target.get_sequence(), str(best_distance), pSNR52_R,
                             scaffold_F, best_target.get_sense(), sense_vs_gene]))
        else:
            print("\t".join([gene] + ['NA'] * 6))

    sys.exit()








