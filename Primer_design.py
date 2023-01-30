#!/usr/bin/env python
#This is python script to design primers to knockout a gene included in file genes.txt from a strain included in file genome.fasta.
#1) It blasts the gene against the genome, outputs the exact allele that the strain has.
#2) extracts 700 base pairs upstream and 700 base pairs downstream of the gene
#3) designs primers to amplify those two fragments using primer3
#4) adds the necessary over hangs for gibson cloning
#5) outputs all of the designed primers with overhangs to a CSV

#Dependencies
#install primer3
#install blast

import os
from Bio import SeqIO
import primer3
import csv

# Create the BLAST database for the genome
os.system("makeblastdb -in genome.fasta -out databaseName -dbtype nucl > /dev/null")

# Perform the BLAST search
os.system("blastn -query genes.txt -db databaseName -outfmt '6' > outfileName")

# Read the BLAST results from the outfile output sequence that has the highest identity and coverage
with open("outfileName", "r") as outfile:
    max_coverage = 0
    max_identity = 0
    max_gene = ""
    for line in outfile:
        fields = line.strip().split("\t")
        qseqid = fields[0]
        sseqid = fields[1]
        identity = float(fields[2])
        coverage = float(fields[3])
        print(max_coverage)
        print(max_identity)
        qstart = int(fields[6])
        qend = int(fields[7])
        sstart = int(fields[8])
        send = int(fields[9])
        print(sseqid)
        if identity >= max_identity and coverage >= max_coverage:
            max_coverage = coverage
            print(max_coverage)
            max_identity = identity
            print(max_identity)
            max_gene = qseqid
            sseqid_record = sseqid
            sstart_record = sstart
            send_record = send


# Extract the aligned gene from the genome
with open("genome.fasta", "r") as genome_file:
    for genome_record in SeqIO.parse(genome_file, "fasta"):
        if genome_record.id == sseqid_record:
            aligned_gene = genome_record.seq[min(sstart_record, send_record) - 1:max(sstart_record, send_record)]
            print(aligned_gene)
# Write the aligned gene to a FASTA file
with open("aligned_gene.fasta", "w") as aligned_gene_file:
    aligned_gene_file.write(">" + max_gene + "\n")
    aligned_gene_file.write(str(aligned_gene))

# Open the genome file as a FASTA file
with open("genome.fasta", "r") as genome_file:
    # Iterate through each record in the FASTA file
    for genome_record in SeqIO.parse(genome_file, "fasta"):
        # Open the gene file as a FASTA file
        with open("aligned_gene.fasta", "r") as gene_file:
            for gene_record in SeqIO.parse(gene_file, "fasta"):
                gene = str(gene_record.seq)
                if gene in str(genome_record.seq):
                    gene_start = str(genome_record.seq).index(gene)
                    gene_end = gene_start + len(gene)
                    if gene_start-700>0:
                        upstream = genome_record.seq[gene_start-700:gene_start-0]
                    else:
                        upstream = genome_record.seq[0:gene_start-0]
                    downstream = genome_record.seq[gene_end:gene_end+700]
                    print(len(upstream))
                    # Write the upstream and downstream sequences to FASTA files
                    with open("upstream.fasta", "w") as upstream_file:
                        upstream_file.write(">upstream_sequence\n")
                        upstream_file.write(str(upstream))
                    with open("downstream.fasta", "w") as downstream_file:
                        downstream_file.write(">downstream_sequence\n")
                        downstream_file.write(str(downstream))
# Open the genome file as a FASTA file
with open("genome.fasta", "r") as genome_file:
    # Iterate through each record in the FASTA file
    for genome_record in SeqIO.parse(genome_file, "fasta"):
        # Open the gene file as a FASTA file
        with open("aligned_gene.fasta", "r") as gene_file:
            for gene_record in SeqIO.parse(gene_file, "fasta"):
                gene = str(gene_record.seq)
                if gene in str(genome_record.seq):
                    gene_start = str(genome_record.seq).index(gene)
                    gene_end = gene_start + len(gene)
                    if gene_start-1000>0:
                        upstream2 = genome_record.seq[gene_start-1000:gene_start-0]
                    else:
                        upstream2 = genome_record.seq[0:gene_start-0]
                    downstream2 = genome_record.seq[gene_end:gene_end+1000]


# Set global parameters for primer3
global_args = {
    'PRIMER_PRODUCT_SIZE_RANGE': [450,650],
    'PRIMER_MIN_SIZE': 25,
    'PRIMER_MAX_SIZE': 25,
    'PRIMER_OPT_SIZE': 25,
    'PRIMER_NUM_RETURN': 1,
    'PRIMER_EXPLAIN_FLAG': '1'
}
global_args2 = {
    'PRIMER_PRODUCT_SIZE_RANGE': [1401,1800],
    'PRIMER_MIN_SIZE': 15,
    'PRIMER_MAX_SIZE': 25,
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_NUM_RETURN': 1,
    'PRIMER_EXPLAIN_FLAG': '1'
}
# Read upstream and downstream sequences from FASTA files
with open('upstream.fasta', 'r') as f:
    upstream_seq = f.read()
with open('downstream.fasta', 'r') as f:
    downstream_seq = f.read()
downstream_seq = downstream_seq.splitlines()
downstream_seq = "\n".join(downstream_seq[1:])
upstream_seq = upstream_seq.splitlines()
upstream_seq = "\n".join(upstream_seq[1:])

#join upstream aligned gene and downstream sequences
complete= str(upstream2+downstream2)
complete2= upstream_seq+aligned_gene+downstream_seq
#downstream_seq = downstream_seq.replace(">downstream_sequence", "")
#upstream_seq = upstream_seq.replace(">upstream_sequence", "")
# Set sequence arguments for primer3
seq_args1 = {
    'SEQUENCE_ID': 'gene1',
    'SEQUENCE_TEMPLATE':  upstream_seq,
    'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [[0, 0, 0, 0]],
    #'SEQUENCE_FORCE_RIGHT_START': (600)
}
seq_args2 = {
    'SEQUENCE_ID': 'gene2',
    'SEQUENCE_TEMPLATE':  downstream_seq,
    'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [[0,0, 0, 0]],
    #'SEQUENCE_FORCE_RIGHT_START': (600)
}
seq_args3 = {
    'SEQUENCE_ID': 'gene3',
    'SEQUENCE_TEMPLATE':  complete,
    'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [[0,0, 0, 0]],
    #'SEQUENCE_FORCE_RIGHT_START': (600)
}
# Run primer3 and get the results
primer_results1 = primer3.bindings.designPrimers(seq_args1, global_args,)
for i in range(primer_results1['PRIMER_PAIR_NUM_RETURNED']):
    left_seq = primer_results1['PRIMER_LEFT_{}_SEQUENCE'.format(i)]

primer_results2 = primer3.bindings.designPrimers(seq_args2, global_args,)
print(primer_results1)
print(len(upstream))
print(primer_results2)

for i in range(primer_results2['PRIMER_PAIR_NUM_RETURNED']):
    right_seq = primer_results2['PRIMER_RIGHT_{}_SEQUENCE'.format(i)]

primer_results3 = primer3.bindings.designPrimers(seq_args3, global_args2,)
print(primer_results3)
for i in range(primer_results3['PRIMER_PAIR_NUM_RETURNED']):
    primer8 = primer_results3['PRIMER_RIGHT_{}_SEQUENCE'.format(i)]
    primer7 = primer_results3['PRIMER_LEFT_{}_SEQUENCE'.format(i)]
#Apramycin sequence
apramycin = "tgtaggctggagctgcttcgaagttcctatactttctagagaataggaacttcggaataggaacttatgagctcagccaatcgactggcgagcggcatcgcattcttcgcatcccgcctctggcggatgcaggaagatcaacggatctcggcccagttgacccagggctgtcgccacaatgtcgcgggagcggatcaaccgagcaaaggcatgaccgactggaccttccttctgaaggctcttctccttgagccacctgtccgccaaggcaaagcgctcacagcagtggtcattctcgagataatcgacgcgtaccaacttgccatcctgaagaatggtgcagtgtctcggcaccccatagggaacctttgccatcaactcggcaagatgcagcgtcgtgttggcatcgtgtcccacgccgaggagaagtacctgcccatcgagttcatggacacgggcgaccgggcttgcaggcgagtgaggtggcaggggcaatggatcagagatgatctgctctgcctgtggccccgctgccgcaaaggcaaatggatgggcgctgcgctttacatttggcaggcgccagaatgtgtcagagacaactccaaggtccggtgtaacgggcgacgtggcaggatcgaacggctcgtcgtccagacctgaccacgagggcatgacgagcgtccctcccggacccagcgcagcacgcagggcctcgatcagtccaagtggcccatcttcgaggggccggacgctacggaaggagctgtggaccagcagcacaccgccgggggtaaccccaaggttgagaagctgaccgatgagctcggcttttcgccattcgtattgcacgacattgcactccaccgctgatgacatcagtcgatcatagcacgatcaacggcactgttgcaaatagtcggtggtgataaacttatcatccccttttgctgatggagctgcattcaaaggccggcattttcagcgtgacatcattctgtgggccgtacgctggtactgcaaatacggcatcagttaccgtgagctgcattttccgctgcataaccctgcttcggggtcattatagcgattttttcggtatatccatcctttttcgcacgatatacaggattttgccaaagggttcgtgtagactttccttggtgtatccaacggcgtcagccgggcaggataggtgaagtaggcccacccgcgagcgggtgttccttcttcactgtcccttattcgcacctggcggtgctcaacgggaatcctgctctgcgaggctggcgggaacttcgaagttcctatactttctagagaataggaacttcgaactgcaggtcgacggatccccggaa"
#Define Primers
primer1 = left_seq
primer2 = upstream_seq[-25:].translate(str.maketrans("ATGC", "TACG"))[::-1]
primer3 = "TGTAGGCTGGAGCTGCTTCGAAGTT"
primer4 = "TTCCGGGGATCCGTCGACCTGCAGT"
primer5 = downstream_seq[:25]
primer6 = right_seq


#Define Overhangs
primer1_overhang = "CCAAGCTTCTCGAGG"
primer2_overhang = "CAGCTCCAGCCTACA"
primer3_overhang = upstream_seq[-15:]
primer4_overhang = downstream_seq[:15].translate(str.maketrans("ATGC", "TACG"))[::-1]
primer5_overhang = "GACGGATCCCCGGAA"
primer6_overhang = "CGGGCTGCAGGAATT"


#Add Overhangs
primer1 = primer1_overhang + primer1
primer3 = primer3_overhang + primer3
primer2 = primer2_overhang + primer2
primer5 = primer5_overhang + primer5
primer4 = primer4_overhang + primer4
primer6 = primer6_overhang + primer6

print(primer1)
print(primer2)
print(primer3)
print(primer4)
print(primer5)
print(primer6)
print(primer7)
print(primer8)
complete=">Complete\n"
complete = str(complete+upstream2+apramycin+downstream2)
print(complete)

primers = [
    {"Primer": "Primer 1", "Sequence": primer1},
    {"Primer": "Primer 2", "Sequence": primer2},
    {"Primer": "Primer 3", "Sequence": primer3},
    {"Primer": "Primer 4", "Sequence": primer4},
    {"Primer": "Primer 5", "Sequence": primer5},
    {"Primer": "Primer 6", "Sequence": primer6},
    {"Primer": "Primer 7", "Sequence": primer7},
    {"Primer": "Primer 8", "Sequence": primer8},
]

with open("primers.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["Primer", "Sequence"])
    writer.writeheader()
    for primer in primers:
        writer.writerow(primer)
