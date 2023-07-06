#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 21:16:41 2022

@author: cmr1
"""

#parses data from SAM formatted alignment output files for graphical presentation using the ggplot2 library

import re
#import pysam
import matplotlib

#filepath for output to be parsed - can be fed from CLI argument too

SAMpath = "/home/cmr1/jhu/project/samples/outXbt2clust.sam"
#the path might have directory information, the filename can be used to identify the source data when naming plots however
filenameonly = re.search("(?<=samples/).+(?=.sam)",SAMpath)
print("processing " + filenameonly.group(0))
#nomatchout = "nomatch_frhit_99.fasta"
accdictpath = "/home/cmr1/jhu/project/samples/noseqs.fasta"

#if not sorted, sort SAM file - this can be disabled by commenting out but could be better as a command line option
#pysam.sort("-o", "plotter_py_sorted_output_temp.sam", SAMpath)

#set up variables for tallying matches - refs[0] is the holder for tallying reads that don't match

nomatches = 0
parse = ""
cols = []
iterator = 0
recorded_something = False
big_meta_list = []

headerref = ""

sam = open(SAMpath,"r")
#nomatch = open(nomatchout,"w")
accdict = open(accdictpath,"r")

#create a dict correlating accession numbers to species names for later
#this assumes fasta format, and works best if all sequence data has been removed previously by the seqstripper.py program
accnames = {}
acc = ""
name = ""
diterator = 0

for header in accdict:
    listed = []
    acc = re.search('(?<=>).+?(?=\s)',str(header))
    name = re.search('(?<=\s).+',str(header))
    if acc == None or name == None:
        continue
    else:
        diterator += 1
        #set up dict to record counted matches, all fragment start positions, and % identity values subdivided by reference genome
        accnames[acc.group(0)] = [name.group(0),0,[],[]]
        
print("dictionary built")

#build a list of lists to collect X and Y data for each of the elements in the dict

never = 0

for line in sam:
    parse = str(line)
    if line[0] == "@":
        continue
    #the format of SAM files is such that all lines beginning with @ will be header lines, and anything below this point will necessarily execute afterward, so I can put everything in one for loop
    else:
        cols=parse.split('\t')
        if never <= 5:
            never = 1
        #base case, increment the no match counter, otherwise search for the right counter to increment
        if cols[2] == '*':
            nomatches+=1
        else:
            #if a match exists to a genome, then that genome will be in the database - if there is no corresponding accession-keyed dict entry, a key error will be thrown
            #increment the match counter
            accnames[cols[2]][1]+=1
            #add a starting coordinate to the list of X coords for the FRP
            accnames[cols[2]][2].append(int(cols[3]))
            #compute percentage identity of the match from CIGAR and MD
            match = (re.findall(r'\d+[A-Z]', cols[5]))
            #counters
            matches = 0
            mismatches = 0
            cigartally = 0
            cigarmiss = 0
            #set a flag just in case the SAM file does not contain an MD tag - in this case we will do our best with CIGAR alone
            mdfound = False
            if match != None:
                for element in match:
                    if element[-1] == "M":
                        cigartally += int(element[:-1])
                    else:
                        cigarmiss += int(element[:-1])
            for column in cols:
                if column[:2] != "MD":
                    continue
                mdfound = True
                match = (re.findall(r'\d+?', column))
                for number in match:
                    matches += int(number)
                    mismatches = len(column) - matches
            #code lifted straight from linked resource "on calculating percent identity from SAM files, except the summed CIGAR M values are used if an MD tag is not available"
            if mdfound == True:
                identity = 100 * (matches / (matches + mismatches))
            else:
                identity = 100 * (cigartally / (cigartally + cigarmiss))
            #error handling
            if identity > 100 or identity < 0:
                identity = -1
            #add it as a y coordinate for the FRP
            accnames[cols[2]][3].append(identity)

print("SAM file parsed, creating plots")

#collect a list of accessions without matches            
notplot = []            
sumofallmatches = 0

for key in accnames:
    #iterate through dictionary and generate plots for every entry with match data
    #skip the ones without any points to plot, collect a list of accession numbers deliberately not plotted to report to the user
    if accnames[key][1] == 0:
        notplot.append(key)
        continue
    #collect a running tally of how many matches there were throughout the data set for later
    sumofallmatches += accnames[key][1]
    #generate labels for the plot and its axes
    print("Plotting accession: " + key)
    matplotlib.pyplot.title("Fragment Recruitment Plot of Matches to " + key)
    plotaccession = str(accnames[key][0])
    matplotlib.pyplot.xlabel("Reference Genome: " + plotaccession)
    matplotlib.pyplot.ylabel("Percent Identity")
    #generate graph axis parameters
    XMin = 0
    #this is probably a waste of resources but it will improve the visual quality of the plots
    XMax = max(accnames[key][2])+1000
    Yfloor = 0
    Ymax = 100
    matplotlib.pyplot.axis([XMin, XMax, Yfloor, Ymax])
    #generate the plotted data
    input1 = matplotlib.pyplot.scatter(accnames[key][2],accnames[key][3],s=1,c="red")
    #display graph in terminal if desired
    matplotlib.pyplot.show()
    plotpath = key + filenameonly.group(0) + ".png"
    matplotlib.pyplot.savefig(plotpath)

print("The following accessions had no matches and were not plotted:\n")
print(notplot)

print("exporting percentages of overall matches to perc_report.txt")
#this is secondary but I wanted to have it because of how overrepresented one strain was in my data set on several occasions"

percpath = str(filenameonly.group(0)) + "perc_report.txt"

perc = open(percpath,"w")
strout = ""

for key in accnames:
    strout = "Accession: " + key + ", " + accnames[key][0] + "\n"
    perc.write(strout)
    strout = str(accnames[key][1]) + " matches out of " + str(sumofallmatches) + "matches found throughout this data set" + "\n"
    perc.write(strout)
    strout = str(100*(accnames[key][1]/sumofallmatches)) + "%\n"
    perc.write(strout)