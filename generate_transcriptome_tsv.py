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
    exons = {}
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            chrom, start, end, strand, name = fields
            enst_id = name.split('_')[0]
            if enst_id not in exons:
                exons[enst_id] = []
            exons[enst_id].append((chrom, int(start), int(end), strand))
    return exons

parser = argparse.ArgumentParser(description="Extract coding sequences from a reference genome using a BED file. NOTE THAT POSITIONS ARE 0 BASED")
parser.add_argument("reference_genome", help="Path to the reference genome in FASTA format.")
parser.add_argument("bed_file", help="Path to the BED file containing coding exon coordinates.")
parser.add_argument("output_tsv", help="Path to the output tsv file containing coding sequences.")

args = parser.parse_args()

print("Parsing reference fasta")

fasta = SeqIO.to_dict(SeqIO.parse(args.reference_genome, "fasta"))

print("Parsing bed file")

exons = parse_bed_file(args.bed_file)

outfile = open(args.output_tsv,"w")

outfile.write("chrom\tpos\tref\talt\tref_trinuc\talt_trinuc\tgene_name\tstrand\tcoding_pos\tcoding_nuc\tcodon_pos\tcodon\talt_codon\tAA_pos\tAA\talt_AA\n")

print("Generating output")

for enst in exons:
    seq = ""
    before_seq = ""
    after_seq = ""
    pos = []
    for exon in exons[enst]:
        chrom, start, end, strand = exon

        seq += fasta[chrom].seq[start:end].upper()
        before_seq += fasta[chrom].seq[(start-1):(end-1)].upper()
        after_seq += fasta[chrom].seq[(start+1):(end+1)].upper()

        pos.extend(list(range(start,end)))
    
    if strand == "+":
        transcript_seq = seq
    else:
        transcript_seq = seq.reverse_complement() 


    for i in range(len(transcript_seq)):
        start = pos[0]
        stop = pos[-1]
        cur_pos = pos[i]
        cur_nuc = seq[i]
        num_nuc = len(transcript_seq)
        num_codon = num_nuc/3
        if strand == "+":
            coding_nuc = transcript_seq[i]
            codon_pos = i % 3
            codon_start = i-codon_pos
            codon_end = i+3-codon_pos
            codon = transcript_seq[codon_start:codon_end]
            AA_pos = math.floor(i/3)
        else:
            coding_nuc = transcript_seq[-(i+1)]
            codon_pos = (num_nuc - i+2) % 3 
            codon_start = -(i+codon_pos+1)
            codon_end = -(i-2+codon_pos) if (i-2+codon_pos) != 0 else None
            codon = transcript_seq[codon_start:codon_end]
            AA_pos = math.floor(num_codon)-math.floor(i/3)

        try:
            AA = codon.translate()
        except BiopythonWarning:
            AA = None
            
        trinuc = before_seq[i] + cur_nuc + after_seq[i]

        for alt_nuc in {"A","C","G","T"} - {cur_nuc}:
            if strand == "-":
                alt_coding_nuc = Seq(alt_nuc).reverse_complement() 
            else:
                alt_coding_nuc = alt_nuc
            alt_codon = codon[:codon_pos] + Seq(alt_coding_nuc) + codon[(codon_pos+1):4]
            try:
                alt_AA = alt_codon.translate()
            except BiopythonWarning:
                AA = None
            alt_trinuc = before_seq[i] + alt_nuc + after_seq[i]
        
            outfile.write(f'{chrom}\t{cur_pos}\t{cur_nuc}\t{alt_nuc}\t{trinuc}\t{alt_trinuc}\t{enst}\t{strand}\t{i}\t{coding_nuc}\t{codon_pos}\t{codon}\t{alt_codon}\t{AA_pos}\t{AA}\t{alt_AA}\n')

#NOTE THAT OUTPUT POSITIONS ARE ZERO BASED!!!
