import subprocess
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio.Blast import NCBIXML

# Define input files
prefix = "barcode24"
query_file = "pFLP.txt"
assembly_file = "/Users/traviskochan/Desktop/Nanopore/assembly_barcode24/00-assembly/draft_assembly.fasta"
output_dir = f"{prefix}_output"  # Use prefix in the output directory name
blast_db = os.path.join(output_dir, f"{prefix}_tmp")  # Use prefix in the database name
output_file = os.path.join(output_dir, f"{prefix}_blast_results.txt")  # Use prefix in the output file name
output_sequence_file = os.path.join(output_dir, f"{prefix}_extracted_sequence.fasta")  # Use prefix in the output sequence file name


# Create BLAST database
makeblastdb_cline = NcbimakeblastdbCommandline(input_file=assembly_file, dbtype="nucl", out=blast_db)
stdout, stderr = makeblastdb_cline()

# Run BLAST command
blast_command = f"blastn -query {query_file} -db {blast_db} -out {output_file} -outfmt 6"
subprocess.run(blast_command, shell=True)

print(f"BLAST results saved to {output_file}")

# Parse BLAST result
with open(output_file, "r") as blast_result:
    lines = blast_result.readlines()
    if not lines:
        print("No BLAST hits found.")
        exit()

    # Extract information about the top hit meeting the criteria
    top_hit = None
    for line in lines:
        fields = line.split("\t")
        alignment_length = int(fields[3])
        if alignment_length >= 980:
            top_hit = fields
            break

    if top_hit is None:
        print("No top hit meeting the criteria found.")
        exit()

    # Extract sequence around the top hit (500bp upstream and 500bp downstream)
    top_hit_start = int(top_hit[8])
    top_hit_end = int(top_hit[9])

    upstream_start = max(0, top_hit_start - 500)
    downstream_end = top_hit_end + 500

    # Extract sequence from the assembly file
    assembly_records = SeqIO.index(assembly_file, "fasta")
    top_hit_sequence = assembly_records[top_hit[1]].seq[upstream_start:downstream_end]

    # Create a SeqRecord object
    top_hit_record = SeqRecord(Seq(str(top_hit_sequence)), id=top_hit[1], description="")

# Write the extracted sequence to a file
    with open(output_sequence_file, "w") as output_sequence:
        SeqIO.write(top_hit_record, output_sequence, "fasta")

# Define input files
query_fasta = output_sequence_file  # Use the output_sequence_file as input
genbank_file = "/Users/traviskochan/Desktop/Assemblies/hvkp5.gb"
output_file = os.path.join(output_dir, f"{prefix}_disrupted_gene.txt")  # Use prefix in the output file name
blast_db = os.path.join(output_dir, f"{prefix}_tmp")  # Use prefix in the database name
fasta_file = "/Users/traviskochan/Desktop/Assemblies/hvkp5.fasta"



# Create BLAST database
makeblastdb_cline = NcbimakeblastdbCommandline(input_file=fasta_file, dbtype="nucl", out=blast_db)
stdout, stderr = makeblastdb_cline()

# Perform local BLAST search
blastn_cline = NcbiblastnCommandline(query=query_fasta, db=blast_db, outfmt=5, out="blast_result.xml")
stdout, stderr = blastn_cline()

# Parse BLAST result and GenBank file
with open("blast_result.xml", "r") as result_handle, open(output_file, "w") as output:
    blast_records = NCBIXML.parse(result_handle)
    genbank_record = SeqIO.read(genbank_file, "genbank")

    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                # Extract information about the hit
                hit_description = alignment.title
                hit_start, hit_end = hsp.sbjct_start, hsp.sbjct_end

                # Find features in the GenBank record that overlap with the hit
                for feature in genbank_record.features:
                    if feature.location.start <= hit_start <= feature.location.end or \
                            feature.location.start <= hit_end <= feature.location.end:
                        # Extract locus tag and gene information from the feature qualifiers
                        locus_tag = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                        gene = feature.qualifiers.get("gene", ["N/A"])[0]

                        # Write information to the output file
                        output.write(f"Hit: {hit_description}\n")
                        output.write(f"Hit Start: {hit_start}\n")
                        output.write(f"Hit End: {hit_end}\n")
                        output.write(f"Locus Tag: {locus_tag}\n")
                        output.write(f"Gene: {gene}\n")
                        output.write("\n")

print(f"Disrupted gene information saved to {output_file}")
