#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 21:16:41 2022

@author: cmr1
"""

#parses data from SAM formatted alignment output files for graphical presentation using the ggplot2 library

import re
import pysam
import matplotlib

#filepath for output to be parsed - can be fed from CLI argument too

SAMpath = "/home/cmr1/jhu/project/samples/fhoutput5.sam"
#nomatchout = "nomatch_frhit_99.fasta"
accdictpath = "/home/cmr1/jhu/project/samples/noseqs.fasta"

#if not sorted, sort SAM file - this can be disabled by commenting out but could be better as a command line option
pysam.sort("-o", "plotter_py_sorted_output_temp.sam", SAMpath)

#set up variables for tallying matches - refs[0] is the holder for tallying reads that don't match

refs = ['No Align']
refcounts = [0]
parse = ""
cols = []
iterator = 0
recorded_something = False
bpx_list = []
idy_list = []
bigmetalist = []

headerref = ""

sam = open("plotter_py_sorted_output_temp.sam","r")
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
        #set up dict to contain a list of all fragment start positions and % identity values subdivided by reference genome
        accnames[acc.group(0)] = [name.group(0),[],[]]
    
print(accnames)

#build a list of lists to collect X and Y data for each of the elements in the dict

for line in sam:
    parse = str(line)
    if line[0] == "@":
        headerref = re.search('(?<=@SQ\tSN:).+(?=\t)',parse)
        if headerref != None:
            refs.append(headerref.group(0))
            refcounts.append(0)        
    #the format of SAM files is such that all lines beginning with @ will be header lines, and anything below this point will necessarily execute afterward, so I can put everything in one for loop
    else:
        cols=parse.split('\t')
        #base case, increment the no match counter, otherwise search for the right counter to increment
        if cols[2] == '*':
            refcounts[0]+=1
            recorded_something = True
            #unidentified sequences stored in separate file for further review
            #nomatch.write(">entry"+str(refcounts[0])+"\n"+cols[9]+"\n")
        else:
            iterator = 0
            recorded_something = False
            #tally up number of matches to assess relative prevalence of a particular reference genome
            for entry in refs:
                if entry == cols[2]:
                    refcounts[iterator]+=1
                    recorded_something = True
                iterator += 1
                if iterator == len(refs) and recorded_something == False:
                    refs.append(cols[2])
                    refcounts.append(1)
                    recorded_something = True
            match = map(str, re.findall(r'\d+[A-Z]', cols[5]))
            if match != None:
                runningmatch = 0
                runningmismatch = 0
                for item in match: 
                    matchcount = 0
                    mismatchcount = 0
                    if item.endswith('M'):
                        matches = int(item[:-1])
                    else:
                        mismatches = int(item[:-1])
                #pull start position of match and add it to the correct list
                fragstart = int(cols[3])
                bpx_list.append(fragstart)
                #compute identity for this match and add it to the correct list
                identity = 100*((matchcount)/float(matchcount + mismatchcount+1))
                idy_list.append(identity)
            
print(refs)
print(refcounts)

print(len(refs)==len(refcounts))

xLabel = "words"
XMin = 0
XMax = 300000
Yfloor = 0
Ymax = 100
AccName = 'Reference Genome: ' + str(xLabel)
matplotlib.pyplot.title('Fragment Recruitment Plot')
matplotlib.pyplot.xlabel(AccName)
matplotlib.pyplot.ylabel('Percent Identity')
matplotlib.pyplot.axis([XMin, XMax, Yfloor, Ymax])

# Generate the scatter plot
input1 = matplotlib.pyplot.scatter(bpx_list,idy_list, c='green', s=0.7, edgecolors='none')

matplotlib.pyplot.show()
outfileName = xLabel + '.png'
matplotlib.pyplot.savefig(outfileName)