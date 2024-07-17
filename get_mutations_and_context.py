import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import warnings
from Bio import BiopythonWarning
import math

warnings.simplefilter("error", BiopythonWarning)

def parse_bed_file(bed_file):
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            chrom, start, end = fields
            regions.append([chrom,int(start),int(end)])
    return regions

parser = argparse.ArgumentParser(description="Extract all possible mutations from a BED file. NOTE THAT POSITIONS ARE 0 BASED")
parser.add_argument("reference_genome", help="Path to the reference genome in FASTA format.")
parser.add_argument("bed_file", help="Path to the BED file containing coordinates.")
parser.add_argument("output_tsv", help="Path to the output tsv file")

args = parser.parse_args()

print("Parsing reference fasta")

fasta = SeqIO.to_dict(SeqIO.parse(args.reference_genome, "fasta"))

print("Parsing bed file")

regions = parse_bed_file(args.bed_file)


outfile = open(args.output_tsv,"w")

outfile.write("chrom\tpos\tref\talt\tref_trinuc\talt_trinuc\n")

print("Generating output")

for region in regions:
    seq = ""
    before_seq = ""
    after_seq = ""
    pos = []
    chrom, start, end = region


    seq = fasta[chrom].seq[start:end].upper()
    before_seq = fasta[chrom].seq[(start-1):(end-1)].upper()
    after_seq = fasta[chrom].seq[(start+1):(end+1)].upper()

    pos.extend(list(range(start,end)))
    

    for i in range(len(seq)):
        start = pos[0]
        stop = pos[-1]
        cur_pos = pos[i]+1 #ONE BASED OUTPUT
        cur_nuc = seq[i]
        num_nuc = len(seq)
            
        trinuc = before_seq[i] + cur_nuc + after_seq[i]

        for alt_nuc in {"A","C","G","T"} - {cur_nuc}:
            alt_trinuc = before_seq[i] + alt_nuc + after_seq[i] 
            outfile.write(f'{chrom}\t{cur_pos}\t{cur_nuc}\t{alt_nuc}\t{trinuc}\t{alt_trinuc}\n')

#NOTE THAT OUTPUT POSITIONS ARE ONE BASED!!!
