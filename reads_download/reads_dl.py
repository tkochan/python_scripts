import subprocess
import os

# Define the folder to save the FASTQ files
output_folder = os.path.expanduser("~/Desktop/Reads")

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

accession_numbers = [
    "SRR25758665","SRR25758666","SRR25758667","SRR25758668","SRR25758669","SRR25758670","SRR25758671","SRR25758672","SRR25758673","SRR25758674","SRR25758675","SRR25758676","SRR25758677","SRR25758678","SRR25758679","SRR25758680","SRR25758681","SRR25758682","SRR25758683","SRR25758684","SRR25758685"
]

# Full path to fastq-dump executable
fastq_dump_path = "/Users/travis.kochan/Desktop/Applications/sratoolkit.3.1.0-mac-x86_64/bin/fastq-dump"

for accession in accession_numbers:
    # Construct the fastq-dump command
    command = f"{fastq_dump_path} --outdir {output_folder} --gzip --readids --read-filter pass --dumpbase --split-3 {accession}"

    # Use subprocess to call the command
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"FASTQ files for {accession} downloaded successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading FASTQ files for {accession}: {e}")
