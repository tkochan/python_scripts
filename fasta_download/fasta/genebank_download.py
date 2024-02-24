import subprocess
import os

# Define the folder to save the FASTA files
output_folder = os.path.expanduser("~/Desktop/Assemblies")

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

accession_numbers = [
    "SAMN35878144", "SAMN35878145", "SAMN35878146", "SAMN35878147", "SAMN35878148",
    "SAMN35878149", "SAMN35878150", "SAMN35878151", "SAMN35878152", "SAMN35878153",
    "SAMN35878154", "SAMN35878155", "SAMN35878156", "SAMN35878157", "SAMN35878158",
    "SAMN11054834", "SAMN35878159", "SAMN35878176", "SAMN35878177", "SAMN35878178",
    "SAMN35878179", "SAMN35878160", "SAMN35878161", "SAMN35878162", "SAMN35878163",
    "SAMN35878164", "SAMN35878165", "SAMN35878166", "SAMN35878167", "SAMN35878168",
    "SAMN35878169", "SAMN35878170", "SAMN35878172", "SAMN35878171", "SAMN35878173",
    "SAMN35878174", "SAMN13161584", "SAMN35878175", "GCF_000755605.1", "SAMD00060934"
]

for accession in accession_numbers:
    # Construct the efetch command
    command = f"esearch -db nucleotide -query {accession} | efetch -format gb"

    # Use subprocess to call the command
    try:
        # Use accession as the filename (remove special characters)
        filename = f"{accession.replace(':', '_').replace('.', '_')}.gb"
        filepath = os.path.join(output_folder, filename)

        with open(filepath, 'w') as file:
            subprocess.run(command, shell=True, check=True, stdout=file)

        print(f"GenBank file for {accession} downloaded successfully. Saved to {filepath}")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading GenBank file for {accession}: {e}")
