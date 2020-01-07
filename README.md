# Candida_CRISPR

This is a series of scripts created to design primers for a CRIPSR knockout screen in the pathogenic fungus Candida albicans.

I designed these to work with a system developed by Min et al. (https://msphere.asm.org/content/1/3/e00130-16, DOI:10.1128/mSphere.00130-16). For specific details, please contact me at andrew.w.pountain@uth.tmc.edu.

In this system, a construct for Cas9 expression is introduced, along with one for gene-targeted guide RNA expression. At the same time, a repair construct containing a selection marker is introduced to cause deletion of the target gene by homologous recombination. Therefore, gene-specific primers are required to make the deletion construct, the guide RNA construct, and validation primers for screening of the deletion.

I've included four scripts:

deletion_primers.py
Takes a fasta file with coding sequences flanked by upstream and downstream regions of the same length (default 1000 bp), along with forward and reverse primer sequences specific for your repair construct. Generates primers with these supplied primer sequences at the 3' end concatenated with 5' sequences to generate the specified length of homology (default 50 bp) of overlap with flanking regions.

sgRNA_primers.py
When supplied with a file containing potential gene-specific guide RNA target sequences, choses the most central of these and uses this to design primers for amplification of a gene-specific guide RNA expression construct.

internal_primers.py
Takes genomic gene sequences and a user-supplied config file, and makes use of the primer3 Python interface (http://primer3.org/manual.html for further details on primer3) to design primer sequences targeting internal sites with a gene. These primers can be used to verify deletion of the desired gene.

make_primer_files.py
Takes the outputs of the above scripts, and combines them into text files for each individual gene, containing all primers required for deletion of that gene.
