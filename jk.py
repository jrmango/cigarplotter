#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 16:20:30 2022

@author: cmr1
"""
from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
import math

# Set initial variables
sampleNumber = 0;
alignedby = 'unknown'
#Create a master list, which becomes a list of lists (one list from each input sample)
master_x_list = []
master_y_list =[]
#master_color_list = []
master_sample_label_list = []

 # Read in each input file one by one
for n in sys.argv[1:]:
    print ("Opening file:", n)
    file_name = n
    sampleNumber += 1
    # At the end of each sample file read, the list of x- and y- coordinates is appended to the corresponding master_list
    x_list = []
    y_list = []
    
    #Example input names: 4e_abau.sam, 4s_abau.sam, is_abau.sam, 4e_abau.psl, 4s_abau.psl, is_abau.psl
    # Gathering information for the legend box
    if "4e" in n:
        legendLabel = "454 even"
    elif "4s" in n:
        legendLabel = "454 staggered"
    elif "is" in n: 
        legendLabel = "Illumina even" 
    else: 
        legendLabel = "unknown"
    
    # Gathering information for the Reference Genome name    
    if "abau" in n:
        xLabel = "Acinetobacter baumannii"
    elif "aodo" in n:
        xLabel = "Actinomyces odontolyticus"
    elif "bcer" in n: 
        xLabel = "Bacillus cereus"
    elif "bvul" in n:
        xLabel = "Bacteroides vulgatus"
    elif "calb" in n:
        xLabel = "Candida albicans"
    elif "ecoli" in n:
        xLabel = "Escherichia coli"
    elif "hpyl" in n:
        xLabel = "Helicobacter pylori"
    elif "msmi" in n:
        xLabel = "Methanobrevibacter smithii"
    elif "sepi" in n:
        xLabel = "Staphylococcus epidermidis"          
    elif "spnu" in n:
        xLabel = "S    #Example input names: 4e_abau.sam, 4s_abau.sam, is_abau.sam, 4e_abau.psl, 4s_abau.psl, is_abau.psl
treptococcus pneumoniae"
    else: 
        xLabel = "unknown"    
        
    # Parse any sam files inputted    
    if file_name.endswith('.sam'):
        print("SAM formatted file found: " + str(file_name))
        alignedby = 'Bowtie2'
        #Parse each file for relevant data    
        with open(file_name) as f:
            for line in f: 
                # Skip header lines
                if line.startswith('@'):
                    continue
                cols = line.rstrip().split()
                if cols[3] == 0:
                    continue
                else:
                    cigar = cols[5]
                    match = map(str, re.findall(r'\d+[A-Z]', cigar))
                    if match:
                        numMatch = 0
                        numMismatch = 0
                        # parse cigar match to find number of matches/mismatches
                        for item in match: 
                            mismatches = 0
                            matches = 0
                            if item.endswith('M'):
                                matches = int(item[:-1])
                                #print ("M item: ", item)
                            else:
                                mismatches = int(item[:-1])
                                #print ("Mismatch item: ", item)
                        
                            numMatch = numMatch + matches
                            numMismatch = numMismatch + mismatches
                            float(numMatch)
                            float(numMismatch)
                            #print ("numMatch:", numMatch, "and mismatch:", numMismatch)
                        
                        # Percent Identity equals number of (matches / (number of matches + mismatches))*100
                        identity = 100*((numMatch)/float(numMatch + numMismatch))    
                        roundX = int(cols[3])
                        #rawX = int(cols[3])
                        #roundX = rawX - (rawX%100)
                        #print (roundX, identity)
                    
                        x_list.append(roundX)
                        y_list.append(identity)
                        
            master_x_list.append(x_list)
            master_y_list.append(y_list)
            master_sample_label_list.append(legendLabel)
                        
    # Parse any psl files inputted
    elif file_name.endswith('psl'):
        print("PSL formatted file found: " + str(file_name))
        alignedby = 'FR-HIT'   
        with open(file_name) as f:
            for line in f:
                cols = line.rstrip().split()
                
                identity = cols[7]
                # Removes the '%' sign
                identity = identity[:-1]
                identity = float(identity)
                
                roundX = cols[9]
                
                x_list.append(roundX)
                y_list.append(identity)
        
                
        master_x_list.append(x_list)
        master_y_list.append(y_list)
    
        master_sample_label_list.append(legendLabel)
    
# Plot Fragment Recruitment Plot    

# This can be altered by the coder, and could be improved to be included as a parameter entered by the user
refLowerRange = 0
refUpperRange = 4000000
lowerIdentity = 80
upperIdentity = 100
# Size of the dots
sz = 1
thisXLabel = 'Reference Genome: ' + str(xLabel)
plt.title('Fragment Recruitment Plot\nAligned via ' + alignedby)
plt.xlabel(thisXLabel)
plt.ylabel('Percent Identity')
plt.axis([refLowerRange, refUpperRange, lowerIdentity, upperIdentity])

# Generate the scatter plot
input1 = plt.scatter(master_x_list[0],master_y_list[0], c='aquamarine', s=sz, edgecolors='none')
input2 = plt.scatter(master_x_list[1],master_y_list[1], c='navy', s=sz, edgecolors='none')
input3 = plt.scatter(master_x_list[2],master_y_list[2], c='magenta', s=sz, edgecolors='none')

# Generate the legend, display it outside of plot. 
x = np.arange(10)
ax = plt.subplot(111)
master_color_list = ['aquamarine', 'navy', 'magenta']
for i in xrange(3):
    ax.plot(x, i * x, color= master_color_list[i], label=master_sample_label_list[i])

# Condense the current height to 85% of original to make room for legend.    
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.85])
# Place the legend                 
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10), fancybox=True, shadow=True, ncol=3 )

plt.show()
outfileName = xLabel + alignedby + '.png'
plt.savefig(outfileName)
