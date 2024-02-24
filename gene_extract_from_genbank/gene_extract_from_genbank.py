from Bio import SeqIO

def extract_gene_sequences(genbank_file, gene_ids, output_dir):
    for gene_id in gene_ids:
        output_file = f"{output_dir}/{gene_id}_sequence.fasta"
        gene_found = False

        # Open the GenBank file and extract the sequence
        with open(genbank_file, "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if "locus_tag" in feature.qualifiers and feature.qualifiers["locus_tag"][0] == gene_id:
                        # Extract the nucleotide sequence
                        gene_sequence = feature.location.extract(record.seq)

                        # Save the sequence to a FASTA file
                        with open(output_file, "w") as output_handle:
                            output_handle.write(f"> {gene_id}\n{gene_sequence}\n")

                        print(f"Nucleotide sequence for gene {gene_id} saved to {output_file}")
                        gene_found = True
                        break

            if not gene_found:
                # If the gene is not found, print a message
                print(f"Gene {gene_id} not found in the GenBank file.")

    # Combine all sequences into a single FASTA file
    combined_output_file = f"{output_dir}/combined_sequences.fasta"
    with open(combined_output_file, "w") as combined_handle:
        for gene_id in gene_ids:
            gene_file = f"{output_dir}/{gene_id}_sequence.fasta"
            try:
                with open(gene_file, "r") as gene_handle:
                    combined_handle.write(gene_handle.read())
            except FileNotFoundError:
                print(f"Skipping {gene_id} for combining as the gene file does not exist.")

    print(f"All sequences combined and saved to {combined_output_file}")

def main():
    # Specify the path to the GenBank file, gene IDs, and output directory
    genbank_file = "/Users/traviskochan/Desktop/Assemblies/hvkp5.gb"
    gene_ids = ["QU725_01075","QU725_17065","QU725_04340","QU725_08205","QU725_25230","QU725_05625","QU725_15125","QU725_20315","QU725_14180","QU725_24470"]
    output_directory = "output_hvKP5"

    # Create the output directory if it doesn't exist
    import os
    os.makedirs(output_directory, exist_ok=True)

    # Extract gene sequences and save to individual and combined FASTA files
    extract_gene_sequences(genbank_file, gene_ids, output_directory)

if __name__ == "__main__":
    main()
