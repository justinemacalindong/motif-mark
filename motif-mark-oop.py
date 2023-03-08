#!/usr/bin/env python

import argparse
import re
import bioinfo
from itertools import product
import cairo

# 1. Parse arguments
def get_args():
    parser = argparse.ArgumentParser(description="Visualize genes and their introns, exons, and motifs in a designated fasta file")
    parser.add_argument("-f", "--file", help="Designates name of fasta file. FILE MUST BE IN SAME FOLDER AS motif-mark-oop.py SCRIPT.", required=True)
    parser.add_argument("-m", "--motifs", help="Designates name of text file of motifs. FILE MUST BE IN SAME FOLDER AS motif-mark-oop.py SCRIPT", required=True)
    return parser.parse_args()

args = get_args()

# Taking name of file to name png and pdf later
file_name = args.file
file_name = file_name.split(".fasta")[0]

# 2. Defining classes
class Gene:
    def __init__(self, start, vert, stop):
        ## Gene class will have start position, stop position, and vertical position for drawing ##
        self.start = start # "0"
        self.stop = stop # length of read
        self.vert = vert # vertical position

    # Methods
    def draw(self, ctx):
        ctx.set_source_rgba(0,0,0) # Set color to black
        ctx.set_line_width(1)
        ctx.move_to(50,self.vert) # x = 50, y = vertical position
        ctx.line_to(self.stop+50,self.vert) # (stop position, vertical position)
        ctx.stroke() # draw line

class Exon:
    def __init__(self, start, vert, stop):
        ## Exon class will have start position, stop position, and vertical position for drawing ##
        self.start = start # starting base position
        self.stop = stop # stop base position
        self.vert = vert # vertical position 

    def draw(self, ctx):
        ctx.set_source_rgba(0,0,0) # Set color to black
        ctx.set_line_width(10)
        ctx.move_to(self.start+50,self.vert) # x = start position, y = vertical position
        ctx.line_to(self.stop+50,self.vert) # x = stop position, y = vertical position
        ctx.stroke() # draw line

class Motif:
    def __init__(self, start, vert, stop):
        ## Motif class will have start position, stop position, and vertical position for drawing ##
        self.start = start # starting base position
        self.stop = stop # stop base position
        self.vert = vert # vertical position

    def draw(self, ctx, colors):
        ctx.set_source_rgba(colors[0],colors[1],colors[2],0.8) # selecting rgb color from list
        ctx.set_line_width(10)
        ctx.move_to(self.start+50,self.vert) # x = start position, y = vertical position
        ctx.line_to(self.stop+50,self.vert) # x = stop position, y = vertical position 
        ctx.stroke() # draw line

def ambig(sequence):
    '''This function takes a motif sequence and replaces nucleotides and ambiguous nucleotides with a regex statement to later search for all motif variations in a read.'''
    ambig_seq = sequence.upper() # capitializes all nucleotides
    ambig_seq = ambig_seq.replace('A', '[Aa]')
    ambig_seq = ambig_seq.replace('T', '[Tt]')
    ambig_seq = ambig_seq.replace('G', '[Gg]')
    ambig_seq = ambig_seq.replace('C', '[Cc]')
    ambig_seq = ambig_seq.replace('U', '[UuTt]')
    ambig_seq = ambig_seq.replace('W', '[AaTtUu]')
    ambig_seq = ambig_seq.replace('S', '[CcGg]')
    ambig_seq = ambig_seq.replace('M', '[AaCc]')
    ambig_seq = ambig_seq.replace('K', '[GgTtUu]')
    ambig_seq = ambig_seq.replace('R', '[AaGg]')
    ambig_seq = ambig_seq.replace('Y', '[CcTtUu]')
    ambig_seq = ambig_seq.replace('B', '[CcGgTt]')
    ambig_seq = ambig_seq.replace('D', '[AaGgTtUu]')
    ambig_seq = ambig_seq.replace('H', '[AaCcTtUu]')
    ambig_seq = ambig_seq.replace('V', '[AaCcGg]')
    ambig_seq = ambig_seq.replace('N', '[AaCcGgTtUu]')
    ambig_seq = ambig_seq.replace('Z', '[-]')
    return ambig_seq # return regex statement


# 4. Clean FASTA and make sequence lines into one line
bioinfo.oneline_fasta(args.file)

# Finding longest read and number of reads in FASTA file to set height width for context
longest_line = 0
line_count = 0 
with open("./one_line.fa", "r") as fasta:
    for line in fasta:
        line_count += 1
        if len(line) > longest_line:
            longest_line = len(line)

width, height = longest_line+100, (line_count/2)*100

# Create context
surface = cairo.PDFSurface(file_name +'.pdf', width, height)
ctx = cairo.Context(surface)
ctx.save()
ctx.set_source_rgb(1,1,1) # setting white background
ctx.paint()

# 5. Parsing through motif file to make usable
orig_motif = [] # list of original motifs
motif_list = [] # list of ambiguous motifs
with open(args.motifs, "r") as motif_file:
    for motif in motif_file:
        motif = motif.strip('\n')
        orig_motif.append(motif)
        motif_list.append(ambig(motif))

# List of colors for drawing motifs (max 5 motifs)
colors = [[0.50, 0.70, 0.52],[0.97, 0.36, 0.01],[0.95,0.77,0.06],[0.42, 0.81, 0.96],[0.99,0.39,0.64]]

# 6. Drawing genes with their exons and motifs one at a time
with open("./one_line.fa", "r") as fasta:

    header_vert = 75 # starting vertical position for header lines
    line_vert = 50 # starting vertical position for gene lines

    for line in fasta:
        line = line.strip() 

        # Header lines
        if line.startswith(">"): # looking for header line
            ctx.move_to(50, header_vert) 
            ctx.set_source_rgba(0,0,0)
            ctx.show_text(line) # showing header line
            header_vert += 75 # move down another 75 pixels

        else:
            # Gene
            gene_line = Gene(0, line_vert, len(line)) # creating gene object with positions
            gene_line.draw(ctx) # draw gene line

            # Exon
            exon_position = re.finditer("[A-Z]+", line)
            for position in exon_position:
                index = position.span() # find start and end positions
            exon_line = Exon(index[0], line_vert, index[1]) # creating exon object with positions
            exon_line.draw(ctx) # draw exon line

            # Motifs
            motif_count = 0 # initializing counter for number of motifs

            # iterating over motif list and drawing all motifs one read at a time
            for motif in motif_list:
                motif_positions = re.finditer(str(motif), line) # find all start and end positions of motif in read

                for position in motif_positions: 
                    color = colors[motif_count] # choose color based on motif count
                    index2 = position.span() # find start and end positions of motif
                    motif_line = Motif(index2[0], line_vert, index2[1]) # creating motif object with position
                    motif_line.draw(ctx,color) # draw motif

                motif_count += 1

            line_vert += 75 # move down 75 pixels

# Drawing legend
legend_count = 0 # initializing counter for colors
start = 0 # initializing counter for start position
stop = 0 # initializing counter for stop position

ctx.set_source_rgba(0,0,0)
ctx.move_to(50+start, height-45)
ctx.show_text("Legend:") # labeling legend

intron_legend = Gene(0, height-30, 12)
intron_legend.draw(ctx) # drawing intron
ctx.move_to(68, height-26)
ctx.show_text("Intron") # labeling intron

exon_legend = Exon(80, height-30, 92)
exon_legend.draw(ctx) # drawing exon
ctx.move_to(148, height-26)
ctx.show_text("Exon") # labeling exon

for motif in orig_motif: # iterating over original motifs
    color = colors[legend_count] # setting motif color
    motif_key = Motif(160+start, height-30, 172+stop)
    motif_key.draw(ctx,color) # drawing motif
    ctx.move_to(228+start, height-26)
    ctx.set_source_rgba(0,0,0)
    ctx.show_text(motif) # labeling motif


    start += 80 # move 80 positions right
    stop += 80 # move 80 positions right
    legend_count +=1 # next color


# Finish drawing
surface.write_to_png(file_name +'.png')
surface.finish()

# Close files
motif_file.close()
fasta.close()