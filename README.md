This is python script to design primers to knockout a gene included in a fastA file from a genome in fastA file.
1) It blasts the gene against the genome, outputs the exact allele that the strain has.
2) extracts 700 base pairs upstream and 700 base pairs downstream of the gene
3) designs primers to amplify those two fragments using primer3
4) adds the necessary over hangs for gibson cloning
5) designs confirmation primers that are designed to be outside of the cloning primers
6) outputs all of the designed primers with overhangs to primers.csv
7) outputs the fragment (upstream-apramycin cassette-downstream) of DNA these primers are meant to amplify and clone as complete.fasta


This requires local installations of primer3 and blast


Usuage:

python3 Primer_design.py {gene fastA file} {genome fastA file}
