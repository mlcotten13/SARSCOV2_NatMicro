#!/usr/local/bin/python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas as pd
import csv
import sys
import os.path

if len(sys.argv)!= 2:
	print("Usage: python map_plot2.py genome_location_table")
	sys.exit()
	
print "Running map_plot2.py really fast"

genome_location_table= sys.argv[1]
outprefix = os.path.splitext(genome_location_table)[0]


fig = plt.figure(figsize=(12, 12))
m = Basemap(projection='gnom', lat_0=1.27, lon_0=32.50, width=740000, height=650000, resolution='h')
m.fillcontinents(color="#e5e4e6", lake_color='#d9eff2')  
m.drawmapboundary(fill_color="#DDEEFF")
m.drawcoastlines()
m.drawcountries()

city_table = csv.reader(open(genome_location_table, "rU"))
city_table.next()#skip header
for row in city_table:	
	x, y = m(float(row[3]),float(row[2]))
	district_label = str(row[0]) +", n="+ str(row[1])
	
	if str(row[4]) == "right":
		district_label = str(row[0]) +"   "
	elif str(row[4]) == "left":
		district_label = "   "+str(row[0])
	plt.text(x-0.1, y+0.1, district_label, fontsize=12, color ='black', horizontalalignment=str(row[4]), 
verticalalignment=str(row[5]))



city_table2 = csv.reader(open(genome_location_table, "rU"))
city_table2.next()
for row in city_table2:
	x, y = m(float(row[3]),float(row[2]))	
	plt.plot(x, y, marker='o', markersize= 15, color = str(row[6]), alpha=0.6)

for pair in [(">10 genomes","#e60033"),("2-10 genomes","#2ca9e1"), ("1 genome", "#4c6473")]:
	plt.scatter([], [], c=str(pair[1]), alpha=1, s=200, label=str(pair[0]))
plt.legend(scatterpoints=1, frameon=False, labelspacing=1, loc='upper left', fontsize=12);

plt.savefig(outprefix+'_map.pdf')
plt.savefig(outprefix+'_map.jpg', dpi =300)

