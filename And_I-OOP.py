#!/usr/bin/env python



import cairo
import math
import re
import random
import numpy as np
import argparse

#Argparse takes in arguments that are callable when running the code.

def args():
    ''' Argparse takes various arguments that will be specific to the users input. Argparse requires
        a fasta file, and motif file. '''
    parser=argparse.ArgumentParser(description = "visualize your motifs!")
    parser.add_argument("-i","--input", help="fasta file with sequences", required = True, type = str)
    parser.add_argument("-m","--motif", help="File containing motifs", required = True, type = str)
    return parser.parse_args()
args = args()

#set argparse 
input_file = args.input
motif_txt = args.motif


''' Dictionaty containing regex expression for IUPAC values'''
IUPAC = {
    "A":"[Aa]",
    "C":"[Cc]",
    "G":"[Gg]",
    "T":"[TtUu]",
    "U":"[UuTt]",
    "W":"[AaTtUu]",
    "S":"[CcGg]",
    "M":"[AaCc]",
    "K":"[GgTtUu]",
    "R":"[AaGg]",
    "Y":"[CcTtUu]",
    "B":"[CcGgTt]",
    "D":"[AaGgTtUu]",
    "H":"[AaCcTtUu]",
    "V":"[AaCcGg]",
    "N":"[AaCcGgTtUu]",
    "Z":"[-]",
}




def parse_motif(motif_txt):
    ''' Extract motifs from motif file '''
    with open(motif_txt, "r") as motif:
        motif_list = [] 
        for line in motif:
            line = line.strip() 
            motif_list.append(line)

    return motif_list



def long_gene():
    ''' Find the longest gene value for reference to plot width  '''
    long_johns = 0
    for i in my_dict:
        val = my_dict[i]
        if len(val) > long_johns:
            long_johns = len(val)

    return long_johns



my_dict = {} #set main dict
''' Retrieving header and sequence information from the fasta '''
with open(input_file, "r") as fa:
    my_dict = {} #fill the fasta dict header key and seq as value
    for line in fa:
        line = line.strip()
        #grab header
        if line[0] == '>':
            header = line
            seq = ''
	#grab seq        
	elif line[0] != '>':
            seq += line
            my_dict[header] = (seq) 


#set the longest gene length 
width_value = long_gene()

#set width and heigth values
width = int(width_value)*1.25
height = int(len(my_dict.values()) * 120) + 20

surface = cairo.SVGSurface("motif_OOP_plot.svg", width, height)
context = cairo.Context(surface)

#set vertical and horizontal positons 
cord_start = 25
position_vert = 25

#extract motifs
motif_list = parse_motif(motif_txt)

def col():
    ''' set values for rgb input'''
    colors = {}
    count = 0
    for i in range(50):
        count += 1
        x,y,z = random.random(),random.random(), random.random()
        colors[count] = [x,y,z]
    return colors

#get color dictionary 
col = col()


class Gene:
    def __init__(self, my_dict):
        self.gene_seq = self.get_gene(my_dict)

    def get_gene(self, my_dict):
        gene_seq = []
        for i in my_dict:
            gene_seq.append(my_dict[i])        
        return gene_seq


class Exon:
    def __init__(self, my_dict):
        self.exon = self.get_exon(my_dict)

    def get_exon(self,my_dict):
        exon_dict = {}
        for i in my_dict.keys():
            exon_seq = my_dict[i] #match exon sequences in main seq dictionary 
            iterate = re.finditer("[A-Z]+", exon_seq) 
            for match in iterate:
                exon_dict[i] = match.span() #finding start and stop place of exon

        return exon_dict


class Motif: 
    def __init__(self,my_dict,motif_list):
        self.motif = self.get_motif(my_dict,motif_list)

    def get_motif(self,my_dict,motif_list):
        allmotes = []
        for i in my_dict:
            gene_seq = my_dict[i]
            mtfL = []
            aladdin = []
            for st in motif_list:
                new_motif = "" #set empty string
                st = st.upper() #convert all to upper
                for char in st:
                    new_motif += str(IUPAC[char]) #add mtf values in IUPAC terms before appendign to list
                mtfL.append(new_motif)
                for j in mtfL:
                    yutter = re.finditer(j,gene_seq) #find mtf in seq
                    monkey = []
                    for matchseq in yutter:
                        mtf_finder = matchseq.span() #gather location of motif
                        monkey.append(mtf_finder)
                aladdin.append(monkey)
            allmotes.append(aladdin)
        return(allmotes)
    


g = Gene(my_dict)
e = Exon(my_dict)
m = Motif(my_dict, motif_list)

motif_count = 0
gene_count = 0

for gene_seq in my_dict.keys(): #set gene_seqs

    #get physical, linear placement of the genes
    context.set_line_width(5)
    context.set_source_rgb(0,0,0)
    context.move_to(cord_start,position_vert +25)
    context.line_to(cord_start+len(g.gene_seq[gene_count]), position_vert + 25) #plot line of gene length
    context.stroke()

    head = list(my_dict.keys())[gene_count]
    context.move_to(cord_start + 10, (position_vert) - 5) #place gene pos
    context.set_source_rgb(0,0,0) #set color to black
    context.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(12)
    context.show_text(head) #show gene name
    context.stroke()
# position_vert += 100
    exon_cord = e.exon[gene_seq]#get exon pos
    exon_cord = list(exon_cord)
    context.rectangle(exon_cord[0], position_vert, exon_cord[1] - exon_cord[0],50) #place exon as large rectangle
    context.set_source_rgb(0.8,0.8,0.8) #set gray scale
    context.fill()

    color_count = 1
    for motif in m.motif[gene_count]:
        for mot in motif:

            x,y,z = col[color_count] #extract values for rgb 

            context.set_line_width(5)
            context.set_source_rgb(x,y,z) #insert rbg values
            context.rectangle(mot[0] + 25, position_vert,mot[1] - mot[0], 50) #place vertical motif rectangle
            context.fill()
        color_count += 1

    position_vert += 100

    gene_count += 1
    motif_count += 1



#get visual key for the motifs
motif_count = 1
lenlg = 415
for motif in motif_list:

    context.set_source_rgb(0,0,0)
    context.move_to(25, lenlg + 5)
    context.show_text(motif) #place motif
    context.set_line_width(1) #draw line between motifs in key
    context.line_to(width - 567, lenlg + 5)
    context.stroke()

    x,y,z = col[motif_count] #extract rgb equivalents in 

    context.set_source_rgb(x,y,z) #set legend colors
    context.rectangle(width - 580, lenlg - 10 , 10, 10) #show motif color
    context.fill()

    motif_count += 1
    lenlg += 22

surface.finish() #fin
