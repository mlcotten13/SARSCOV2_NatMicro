#!/usr/local/bin/python
import sys
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
import sys
import os.path 
from Bio.Blast import NCBIXML
from Bio import Entrez
import time
import csv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from collections import defaultdict
 
if len(sys.argv)!= 6:
	print("Usage: python HiLiter_aa2.py sequences.fasta feature_table.csv annotation_table.csv protein annotate_cutoff")
	sys.exit()
	
print "Running HiLiter_aa2.py really fast"
all_seqs= sys.argv[1]
outprefix = os.path.splitext(all_seqs)[0]
feature_table = sys.argv[2]
annotation_table = sys.argv[3]
protein_name = sys.argv[4]
annotate_cutoff = int(sys.argv[5])
seq_names =[]
sequences =[]

for record in SeqIO.parse(open(all_seqs, "rU"), "fasta"):
	seq_names.append(record.id)
	sequence_string = str(record.seq)	
	sequences.append(sequence_string)
total_number_sequences = len(sequences)
genome_length = len(sequences[0])
print "total_number_sequences" 
print total_number_sequences

if total_number_sequences >0 and total_number_sequences <=30:
	ytick_label_size = 5
elif total_number_sequences >30 and total_number_sequences <=100:
	ytick_label_size = 4
else:
	ytick_label_size = 1

genome_length = len(sequences[0])
print "genome_length" 
print genome_length

if genome_length >= 0 and genome_length < 200:
	marker_size = 0.5
elif genome_length >= 200 and genome_length < 350:
	marker_size = 3	
elif genome_length >= 350 and genome_length < 600:
	marker_size = 6
elif genome_length >= 600 and genome_length < 1500:
	marker_size = 10
elif genome_length >= 1200 and genome_length < 3000:
	marker_size = 12
else:
	marker_size = 15
print "marker_size" 
print marker_size
	
graph_pad = 5 * marker_size

graph_names=[]
graph_names.append(seq_names[0])

def vergleichen(seqA, seqB, nameB, index):	
	for i in range (len(seqA)):
		if seqB[i]!=seqA[i] and seqB[i] == "D":
			color = "orangered"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "E":
			color = "crimson"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "K":
			color = "teal"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "R":
			color = "#2f5d50"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "H":
			color = "#38b48b"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "Q":
			color = "slateblue"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "N":
			color = "seagreen"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "P":
			color = "indigo"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "F":
			color = "darkblue"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "W":
			color = "coral"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "G":
			color = "#006e54" #green onion color
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "A":
			color = "tan"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "L":
			color = "#59b9c6"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "I":
			color = "lightseagreen"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "V":
			color = "brown"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "S":
			color = "indianred"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "T":
			color = "darksalmon"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "Y":
			color = "tomato"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "C":
			color = "cadetblue"
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqB[i]!=seqA[i] and seqB[i] == "M":
			color = "skyblue"
			line = index
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			diff_list.append((i+1,line,color))	
			complete_change_list.append(i+1)
		elif seqB[i]=="-":
			color = "#c0c6c9"
			line = index
			diff_list.append((i+1,line,color))	
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
		else:
			continue
diff_list=[]
table_list=[]
complete_change_list = []
 
for i in range(1, len(sequences)):
	vergleichen(sequences[0], sequences[i], seq_names[i], i)

all_changes=[]
difference_summary = open(outprefix +"_all_changes.csv", "a")
for element in table_list:
	row_number =  str(element[0])
	sample = element[1]
	position = str(element[2])
	change_annotation = str(element[3]+position+element[4])
	print_string = str(sample+","+row_number+","+position+","+change_annotation)
	print >> difference_summary, print_string
	all_changes.append(change_annotation)
difference_summary.close()

occurrence = defaultdict(lambda: 0)

for element in all_changes:
	occurrence[element] += 1

variant_catalog = open(outprefix +"_variant_catalog.csv", "w")
print_string_header = "Position,Variant,From,To,Count"
print >> variant_catalog, print_string_header
variant_catalog.close()
for key, value in occurrence.items():
	position = str(key)[1:-1]
	print_string = position+","+str(key)+","+str(key[0])+","+str(key[-1])+","+str(value)
	variant_catalog = open(outprefix +"_variant_catalog.csv", "a")
	print >> variant_catalog, print_string
	variant_catalog.close()

annotation_list=[]
frequency_count = defaultdict(int)
for item in all_changes:
	frequency_count[item] += 1
for key, value in frequency_count.items():
	if "-" in key:
		pass
	elif value >= annotate_cutoff:
		print key, value	
		annotation_list.append(key)

def getKey(item):
	return item[0]
	
simplified_annotation_table = open(outprefix +"_annotations.csv", "a")
for annotation in annotation_list:
	sorting_list = []
	for element in table_list:
		if annotation in element:
			sorting_list.append(element)
	ordered_sorting_list= sorted(sorting_list, key=getKey)
	print ordered_sorting_list 
	print ordered_sorting_list[-1]						
	print_string = ', '.join(str(elem) for elem in ordered_sorting_list[-1])
	print >> simplified_annotation_table, print_string	

simplified_annotation_table.close()
ticks=[] 
for x in range(len(sequences)+1):
	ticks.append(x)
fig = plt.figure()

ax2 = plt.subplot2grid((28,1), (6,0), rowspan = 22)
ax2.set_ylim(0.9, len(sequences)-1) 
ax2.set_xlim(-graph_pad,(genome_length+graph_pad))
ax2.set_yticks(ticks) #to add sequenc elabels to graph
ax2.set_yticklabels(seq_names, size = ytick_label_size) #to add sequence labels to graph
ax2.grid(False)
ax2.set_xlabel(protein_name+' Protein Position')

for i in range(len(diff_list)):
	color = diff_list[i][2]
	ax2.broken_barh([(int(diff_list[i][0]), marker_size)] , (int(diff_list[i][1]), 0.85), facecolors=(color),edgecolors=('None')) 

annotation_table_rdr = csv.reader(open(annotation_table, "rU"))
for row in annotation_table_rdr:
	x_pos = int(row[2])+15
	y_pos = float(row[0])
	change = str(row[5])
	ax2.annotate(change, xy=(1,1), xytext=(x_pos, y_pos), fontsize=6)	

feature_table = csv.reader(open(feature_table, "rU"))
feature_names =[]
feature_map_positions = []
feature_lengths = []
feature_colors=[]
feature_row = []
feature_table.next()
for row in feature_table:	
	feature_name = str(row[0])
	feature_names.append(feature_name)
	feature_start = int(row[1])
	feature_map_positions.append(feature_start)
	feature_length = int(row[2])- int(row[1])
	feature_lengths.append(feature_length)
	feature_color = str(row[4])
	feature_colors.append(feature_color)
	feature_row.append(int(row[3]))
	
ax3 = plt.subplot2grid((28,1), (1,0), rowspan = 5)
ax3.set_ylim(0.8, 6.2) 
ax3.set_xlim(-graph_pad,(genome_length+graph_pad))
ax3.set_xticklabels([])
ax3.set_yticklabels([])

ax3.grid(False)
ax3.get_yaxis().set_visible(False)
for i in range(len(feature_map_positions)):
	ax3.broken_barh([(int(feature_map_positions[i]), int(feature_lengths[i]))],(float(feature_row[i]), 0.9), facecolors=(str(feature_colors[i])), edgecolors=('None'))
	ax3.annotate(feature_names[i], xy=(1,1), xytext=(feature_map_positions[i]+int(feature_lengths[i])+2, (float(feature_row[i]))+0.2), fontsize=5)

plt.tight_layout()
plt.savefig(outprefix+'_differences.pdf', dpi=300)
plt.savefig(outprefix+'_differences.jpg', dpi=300)

print "That's All Folks!"
