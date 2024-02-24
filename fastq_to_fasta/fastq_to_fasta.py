from Bio import SeqIO
import gzip

def convert_fastq_to_fasta(input_fastq, output_fasta):
    with gzip.open(input_fastq, "rt") as fastq_file, open(output_fasta, "w") as fasta_file:
        for record in SeqIO.parse(fastq_file, "fastq"):
            SeqIO.write(record, fasta_file, "fasta")

# Replace "barcode01.fastq.gz" with the actual file name
input_fastq_file = "barcode17.fastq.gz"

# Replace "output.fasta" with the desired output file name
output_fasta_file = "output_17.fasta"

convert_fastq_to_fasta(input_fastq_file, output_fasta_file)

print(f"Conversion complete. Fasta file saved as {output_fasta_file}")
