<b>Circular Plasmid Alignment Visualization</b>

This Python script runs BLAST to align multiple plasmid sequences against a reference sequence and visualizes the alignments as a circular map using Matplotlib. It is useful for comparing structural similarities between plasmids and visualizing the presence of specific gene clusters or sequences relative to a reference.

<b>Features</b>

✔ Runs BLAST locally to generate pairwise alignments between a reference and multiple plasmids.
✔ Produces a circular alignment plot showing the reference as a ring and each plasmid as concentric arcs.
✔ Customizable colors, labels, and tick marks for easy interpretation.
✔ Saves the plot as high-resolution PNG and PDF.

<b>Requirements</b>

Python 3.8+

BLAST+ command-line tools
 installed and accessible in your PATH

Python packages:

pip install biopython matplotlib numpy

<b>Usage</b>

Edit the input file paths at the top of the script:

``` reference_file = "directory/to/reference.fasta"
plasmid_files = [
    "directory/to/plasmid1.fasta",
    "directory/to/plasmid2.fasta",
    "directory/to/plasmid3.fasta"
]
blast_result_file = "directory/to/alignments.txt"
```

Run the script:

``` python circular_alignment.py ```


<b>Output:</b>

BLAST alignment results saved as alignments.txt

Circular plot saved as:

circular_alignment2.png

circular_alignment2.pdf

<b>Output Visualization</b>

Reference sequence is displayed as a full circle in light blue.

Each plasmid is shown as a colored arc aligned to the reference.

Tick marks indicate positions in kilobases.

A legend identifies plasmids by color.

The reference name and size are displayed in the center of the plot.

<b>Example</b>

If your reference is pK2044 (200 kb) and you align 3 plasmids, the resulting figure will show:

Outer arcs representing the plasmids aligned to their corresponding regions on the reference.

Tick labels in kb.

Legend on the right with plasmid names.

<b>Notes</b>

Make sure BLAST+ is installed and available in your system path. You can test by running:

``` blastn -version ```


For large genomes or many plasmids, increase available memory or adjust BLAST parameters.

Customize colors and plot settings inside the script as needed.
