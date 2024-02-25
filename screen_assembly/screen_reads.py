from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio.Blast import NCBIXML


# Define the input files
fasta_file = "/Users/traviskochan/Desktop/Nanopore/output_17.fasta"
blast_output_file = "output_17.txt"
output_file = "pSamhyg.fasta"
min_alignment_length = 1000

# Read blast output to get contig names and their alignment lengths
identified_contigs = {}
with open(blast_output_file, "r") as blast_file:
    for line in blast_file:
        fields = line.split("\t")
        contig_name = fields[1]
        alignment_length = int(fields[3])
        if alignment_length >= min_alignment_length:
            identified_contigs[contig_name] = alignment_length

# Filter and write sequences to the output file
with open(output_file, "w") as output:
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in identified_contigs:
            SeqIO.write(record, output, "fasta")

print(f"Sequences saved to {output_file}")

# Define input files
query_fasta = "pSamhyg.fasta"
genbank_file = "/Users/traviskochan/Desktop/Assemblies/hvkp5.gb"
output_file = "disrupted_gene.txt"
blast_db = "tmp"
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
