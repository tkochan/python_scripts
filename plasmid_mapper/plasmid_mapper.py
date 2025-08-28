import subprocess
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np

# Top section: specify your reference and plasmid file paths
reference_file = "directory/to/reference.fasta"  # Path to your reference FASTA file
plasmid_files = ["directory/to/plasmid1.fasta",  # List multiple plasmid FASTA files here
                 "directory/to/plasmid2.fasta",
                 "directory/to/plasmid3.fasta"
                 
                 ]
blast_result_file = "/Users/travis.kochan/Desktop/alignments.txt"  # Path to your BLAST alignment result

# Function to run BLAST
def run_blast(reference, plasmids, output_file):
    """Runs BLAST and saves the output in tabular format."""
    db_name = "blast_db"

    # Create BLAST database from reference
    subprocess.run(["makeblastdb", "-in", reference, "-dbtype", "nucl", "-out", db_name], check=True)

    # Run BLAST for each plasmid
    with open(output_file, "w") as out:
        for plasmid in plasmids:
            subprocess.run([
                "blastn", "-query", plasmid, "-db", db_name,
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
            ], stdout=out, check=True)

# Step 1: Function to read sequences from a FASTA file
def read_fasta(file):
    """Reads sequences from a FASTA file."""
    with open(file, "r") as f:
        return list(SeqIO.parse(f, "fasta"))

# Step 2: Function to parse BLAST output (tabular format)
def parse_blast_result(blast_file):
    """Parse BLAST output in tabular format."""
    alignments = []
    with open(blast_file, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            plasmid_id, ref_id, _, _, _, _, q_start, q_end, s_start, s_end, _, _ = fields[:12]
            alignments.append((plasmid_id, ref_id, int(q_start), int(q_end), int(s_start), int(s_end)))
    return alignments

# Step 3: Function to create a circular plot of the plasmid alignments to the reference
def plot_circular_alignment(reference, plasmids, alignments):
    """Generates a circular plot of the plasmid alignments to the reference."""

    # Create a figure for the circular plot
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    ax.set_theta_offset(np.pi / 2)  # Set zero at the top
    ax.set_theta_direction(-1)      # Reverse the direction (clockwise)

    # Remove radial grid lines and labels
    ax.set_yticklabels([])  # Remove radial labels
    ax.grid(False)  # Remove the grid lines

    # Plot the reference as a full circle at a smaller radius to decrease space
    ref_len = len(reference)
    ref_radius = 1  # Set reference plot radius
    ref_linewidth = 7  # Define reference line thickness (to match plasmids)
    ax.plot(np.linspace(0, 2*np.pi, ref_len), [ref_radius]*ref_len, linewidth=ref_linewidth, label="pK2044", color="lightblue")  # Light blue color for reference

    # Plot each plasmid's alignment at different radii (incrementing for each plasmid)
    alignment_radius_increment = 0.06  # Increment radius for each plasmid
    unique_labels = set()

    # Define custom colors: red, purple, and blue
    custom_colors = ["red", "black", "grey", "brown", "orange", "blue", "green", "pink", "purple"]

    # Create a dictionary to group alignments by plasmid
    plasmid_alignments = {}
    for plasmid, _, q_start, q_end, s_start, s_end in alignments:
        if plasmid not in plasmid_alignments:
            plasmid_alignments[plasmid] = []
        plasmid_alignments[plasmid].append((q_start, q_end, s_start, s_end))

    # Plot each plasmid's alignment with custom colors
    alignment_radius = 1.06  # Start the first plasmid slightly above the reference circle
    for i, (plasmid, align_list) in enumerate(plasmid_alignments.items()):
        color = custom_colors[i % len(custom_colors)]  # Assign the custom colors in a cycle
        label_added = False  # Flag to check if the label is added for the first time
        for q_start, q_end, s_start, s_end in align_list:
            s_start, s_end = min(s_start, s_end), max(s_start, s_end)  # Ensure positive range

            # Increase number of points to make the arc smoother
            num_points = 500  # Increase this for a smoother appearance
            theta = np.linspace(0, 2*np.pi, num_points)

            # Ensure s_end does not exceed the maximum index
            s_end = min(s_end, ref_len - 1)

            # Plot the alignment arc for each plasmid
            region_start = int(np.linspace(0, num_points, num_points)[int(s_start / ref_len * num_points)])
            region_end = int(np.linspace(0, num_points, num_points)[int(s_end / ref_len * num_points)])
            ax.plot(theta[region_start:region_end], [alignment_radius] * (region_end - region_start), linewidth=7, color=color)

            # Add the legend label only once for each plasmid
            if not label_added:
                ax.plot([], [], label=plasmid, color=color, linewidth=7)  # Invisible plot to add label only once
                label_added = True

        alignment_radius += alignment_radius_increment  # Increment the radius for the next plasmid

    # Adjust limits and add title
    ax.set_ylim(0, alignment_radius + 0.05)  # Adjust the upper limit dynamically based on the last plasmid
    ax.set_title("", va='bottom')

    # Set the ticks and labels to represent kilobases (kb)
    tick_positions = np.linspace(0, 2*np.pi, ref_len, endpoint=False)
    tick_labels = [f"{int(i / 1000)} kb" for i in range(0, ref_len, int(ref_len / 10))]  # Assuming reference is in base pairs

    # Remove the last tick (0kb) to prevent overlap
    ax.set_xticks(tick_positions[::int(ref_len / 10)][: -1])  # Exclude the last tick
    ax.set_xticklabels(tick_labels[:-1])  # Exclude the last label

    # Add short tick lines
    tick_line_length = 0.05  # Adjust length of the ticks (make them shorter)
    for tick in tick_positions[::int(ref_len / 10)][: -1]:
        ax.plot([tick, tick], [ref_radius, ref_radius + tick_line_length], color='black', linewidth=0.7)  # Short line at each tick

    # Move the tick labels inside the reference ring
    for label in ax.get_xticklabels():
        label.set_fontsize(10)  # Set a smaller font size for the labels
        label.set_verticalalignment('center')  # Center-align the labels vertically
        label.set_horizontalalignment('center')  # Center-align the labels horizontally
        label.set_position((label.get_position()[0], label.get_position()[1] + 0.55))  # Move the labels inside the reference ring

    # Hide the outer border (black ring)
    ax.spines['polar'].set_visible(False)  # Hide the polar spine (outer border)

    # Add the reference plasmid name and size in the middle of the ring
    reference_plasmid_name = "pK2044"  # Example reference name
    reference_plasmid_size_kb = ref_len / 1000  # Size in kilobases
    ax.text(0, 0, f"{reference_plasmid_name}\n{int(reference_plasmid_size_kb)}kb",
            horizontalalignment='center', verticalalignment='center', fontsize=14) #, fontweight='bold')

    # Move legend off to the side
    # Move legend closer to the circle
    plt.legend(bbox_to_anchor=(1.02, 0.8), loc='center left', frameon=False)
    # Save the figure as a high-resolution PNG before showing it
    plt.savefig("circular_alignment2.png", dpi=300, bbox_inches='tight')
    plt.savefig("circular_alignment2.pdf", bbox_inches='tight')
    plt.show()


# Step 4: Main function to read files and generate the plot
def main():
    # Run BLAST alignment
    run_blast(reference_file, plasmid_files, blast_result_file)
    # Step 4a: Read the reference and plasmid sequences
    reference = read_fasta(reference_file)[0]
    plasmids = [read_fasta(plasmid_file) for plasmid_file in plasmid_files]

    # Step 4b: Parse the BLAST alignment result
    alignments = parse_blast_result(blast_result_file)

    # Step 4c: Generate and plot the circular alignment
    plot_circular_alignment(reference, plasmids, alignments)

# Run the script
if __name__ == "__main__":
    main()
