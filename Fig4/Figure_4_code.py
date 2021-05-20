#!/usr/local/bin/python
import sys
import csv
import sys
import os
from pylab import *
import numpy as np
import numpy.random
import matplotlib.pyplot as plt
from pylab import imshow, show, get_cmap
import matplotlib.cm as cm
import pylab as pl
import matplotlib as mpl

print('Running variance_barplot1.py')

if len(sys.argv)!=4:
	print("Usage: python variance_barplot1.py input_table.csv outprefix_short y_scale_max")
	sys.exit()

input_table = sys.argv[1]
outprefix = os.path.splitext(input_table)[0]
outprefix_short = sys.argv[2]
y_scale_max = float(sys.argv[3])

color_dictionary = {"ORF8":"#00035b","nsp6TM":"#007bbb","ORF9N":"#0485d1","nsp15endo":"#38b48b","ORF7a":"#00a497","nsp1":"#485690","RDRP":"#4d5aaf","nsp9RNAbp":"#59b9c6","nsp4":"#80aba9","nsp7":"#5c9291","ORF4E":"#6c848d","nsp16OMT":"#95949a","nsp2":"#a2d7dd","nsp14ExoN":"#a8bf93","nsp3":"#b00149","ORF3a":"#c8c2be","nsp10CysHis":"#e2041b","nsp13hel":"#ec6d71","ORF5M":"#ed6d3d","nsp53Cpro":"#ee7948","nsp8Rep":"#fcd575","Spike":"#fd3c06","ORF6": "#4f455c"}


index_list =[]
variance_list = []
mean_list =[]
color_list =[]
legend_list =[]
one_minus_mean_list = []
starting_csv = open(input_table, "rU")
reader = csv.reader(starting_csv)
next(reader)
for row in reader:	
	index_list.append(int(row[0]))
	variance_list.append(float(row[2]))	
	mean_list.append(float(row[3]))
	one_minus_mean_list.append(float(row[4]))
	this_color = color_dictionary.get(str(row[6])) 	
	color_list.append(this_color)
	protein_color_position=str(row[6])+","+str(this_color)+","+str(row[5])
	legend_list.append(protein_color_position)

legend_list_unique = list(set(legend_list))
print(legend_list_unique) 

fig = plt.figure(figsize=(15,4))
ax2 = fig.add_subplot(211)
ax2.set_ylim(0, y_scale_max)
ax2.set_xlim(0,385)
ax2.set_xticklabels([])
title(outprefix_short+"_domain_1-bitscore_mean")
plt.grid(True)
plt.xlabel('Genome position')
plt.ylabel('1-bitscore_mean')

plt.bar(index_list, one_minus_mean_list, color = color_list, width=1, alpha = 0.8)

p0 = Rectangle((0, 0), 1, 1, fc="#485690")
p1 = Rectangle((0, 0), 1, 1, fc="#a2d7dd")
p2 = Rectangle((0, 0), 1, 1, fc="#b00149")
p3 = Rectangle((0, 0), 1, 1, fc="#80aba9")
p4 = Rectangle((0, 0), 1, 1, fc="#ee7948")
p5 = Rectangle((0, 0), 1, 1, fc="#007bbb")
p6 = Rectangle((0, 0), 1, 1, fc="#5c9291")
p7 = Rectangle((0, 0), 1, 1, fc="#fcd575")
p8 = Rectangle((0, 0), 1, 1, fc="#59b9c6")
p9 = Rectangle((0, 0), 1, 1, fc="#e2041b")
p10 = Rectangle((0, 0), 1, 1, fc="#4d5aaf")
p11 = Rectangle((0, 0), 1, 1, fc="#ec6d71")
p12 = Rectangle((0, 0), 1, 1, fc="#a8bf93")
p13 = Rectangle((0, 0), 1, 1, fc="#38b48b")
p14 = Rectangle((0, 0), 1, 1, fc="#95949a")
p15 = Rectangle((0, 0), 1, 1, fc="#fd3c06")
p16 = Rectangle((0, 0), 1, 1, fc="#c8c2be")
p17 = Rectangle((0, 0), 1, 1, fc="#6c848d")
p18 = Rectangle((0, 0), 1, 1, fc="#ed6d3d")
p19 = Rectangle((0, 0), 1, 1, fc="#4f455c")
p20 = Rectangle((0, 0), 1, 1, fc="#00a497")
p21 = Rectangle((0, 0), 1, 1, fc="#00035b")
p22 = Rectangle((0, 0), 1, 1, fc="#0485d1")

plt.legend([p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22], ["nsp1","nsp2","nsp3","nsp4","nsp53Cpro","nsp6TM","nsp7","nsp8Rep","nsp9RNAbp","nsp10CysHis","RDRP","nsp13hel","nsp14ExoN","nsp15endo","nsp16OMT","Spike","ORF3a","ORF4E","ORF5M","ORF6","ORF7a","ORF8","ORF9N"],loc='upper left', bbox_to_anchor=(-0.01, -0.2), prop={'size':8}, ncol=12)

plt.savefig(outprefix_short+'_1-bitscore_mean.pdf')
plt.savefig(outprefix_short+'_1-bitscore_mean.jpg', dpi = 300)

print("That's all Folks!")
