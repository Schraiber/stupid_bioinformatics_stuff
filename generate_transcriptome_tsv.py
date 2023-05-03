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

parser = argparse.ArgumentParser(description="Extract coding sequences from a reference genome using a BED file.")
parser.add_argument("reference_genome", help="Path to the reference genome in FASTA format.")
parser.add_argument("bed_file", help="Path to the BED file containing coding exon coordinates.")
parser.add_argument("output_tsv", help="Path to the output tsv file containing coding sequences.")

args = parser.parse_args()

print("Parsing reference fasta")

fasta = SeqIO.to_dict(SeqIO.parse(args.reference_genome, "fasta"))

print("Parsing bed file")

exons = parse_bed_file(args.bed_file)

outfile = open(args.output_tsv,"w")

outfile.write("chrom\tpos\tref\tref_trinuc\tgene_name\tstrand\tcoding_pos\tcoding_nuc\tcodon_pos\tcodon\tAA_pos\tAA\n")

print("Generating output")

for enst in exons:
    seq = ""
    before_seq = ""
    after_seq = ""
    pos = []
    for exon in exons[enst]:
        chrom, start, end, strand = exon

        seq += fasta[chrom].seq[start:end]
        before_seq += fasta[chrom].seq[(start-1):(end-1)]
        after_seq += fasta[chrom].seq[(start+1):(end+1)]

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
            print(enst)
            print(strand)
            print(len(transcript_seq))
            print(i)
            print(codon_pos)
            print(codon)
            input()
            AA = None
            
        trinuc = before_seq[i] + cur_nuc + after_seq[i]
        outfile.write(f'{chrom}\t{cur_pos}\t{cur_nuc}\t{trinuc}\t{enst}\t{strand}\t{i}\t{coding_nuc}\t{codon_pos}\t{codon}\t{AA_pos}\t{AA}\n')
        #if (AA == "*") and ((codon_pos == 2 and strand=="+") or (codon_pos == 0 and strand=="-")) : break 
        

