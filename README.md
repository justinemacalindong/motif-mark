# motif-mark

This script, **motif-mark-oop.py**, takes an input FASTA and text file containing motif sequences and outputs a visual PNG and PDF representation of the FASTA file. It labels the gene sequences contained in the FASTA file with their exons and the given motifs. The length of all lines and boxes are to scale and correspond to the length of each sequence. 

## Input files:
- **-f**: Fasta file containing gene sequences
    - Introns - lowercase
    - Exons - uppercase
- **-m**: Text file containing motif sequences

## Output files:
- PNG file named by input FASTA file name
- PDF file named by input FASTA file name

