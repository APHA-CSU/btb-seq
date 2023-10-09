#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#import the required modules
import pandas as pd
import sys
import os

#get the required variables from when python script is run
args=sys.argv

mask = sys.argv[1]
positions =sys.argv[2]
output =sys.argv[3]

#open and read the SNP positions for the sample, get them in a list
file = open(positions).read().split('\n')

file = file[:-1]
pos = []

for x in file:
    pos.append(int(x))

#test if the samples has too many SNPs, if it does it is likey contaminated and will take too long to go through the 10 base SNP filter so output the old mask, will likey be caught by the other filters
"""
if len(pos) > 5000:
    bedfile = pd.DataFrame(columns=["Sample","Start","End"])
    
    f = open(mask).read().split('\n')
    f = f[:-1]
    datadriven = []
    for x in f:
        datadriven.append(x.split("\t"))

    for x in datadriven:
        bedfile.loc[-1] = ["LT708304-Mycobacteriumbovis-AF2122-97",int(x[1]),int(x[2])]
        bedfile.index = bedfile.index + 1
        bedfile = bedfile.sort_index()

    bedfile = bedfile.sort_values("Start")

    bedfile.to_csv(output, sep='\t', index=False, header=None)
    sys.exit()
"""

#calculate if the values are within distance value of the numbers either side of them on the list, for the first and last value on the list, it checks on the value above and below it respectively.
distance = 10
remove = []
counter = 0
for x in pos:
    if counter == 0:
        x_upper = x + distance
        if pos[(counter+1)] < x_upper:
            remove.append(x)
    elif counter == (len(pos) - 1):
        x_lower = x - distance
        if pos[(counter-1)] > x_lower:
            remove.append(x)
    else:
        x_upper = x + distance
        x_lower = x - distance
        if pos[(counter+1)] < x_upper or pos[(counter-1)] > x_lower:
            remove.append(x)
    counter = counter + 1

#Looking to see if the next value in the list is within distance of the previous value, if it is then it will move on, if not, it will end the chain and you will get a start and end value for the chain
counter = 0
chain = []
for x in remove:
    if counter == 0:
        start = x
    x_upper = x + distance
    if counter == (len(remove) - 1):
        end = x
        chain.append([start,end])
        break
    if remove[(counter + 1)] > x_upper:
        end = x
        chain.append([start,end])
        start = remove[(counter + 1)]
    counter = counter + 1

#convert the chain value into something similar to a bedfile
bedfile = pd.DataFrame(columns=["Sample","Start","End"])
for x in chain:
    bedfile.loc[-1] = ["LT708304-Mycobacteriumbovis-AF2122-97",(x[0]-1),x[1]]
    bedfile.index = bedfile.index + 1
    bedfile = bedfile.sort_index()

#import the current mask bedfile
file = open(mask).read().split('\n')
file = file[:-1]
datadriven = []
for x in file:
    datadriven.append(x.split("\t"))

#add the current mask positions to the base SNP filter positions
for x in datadriven:
    bedfile.loc[-1] = ["LT708304-Mycobacteriumbovis-AF2122-97",int(x[1]),int(x[2])]
    bedfile.index = bedfile.index + 1
    bedfile = bedfile.sort_index()

#sort these values from low to high
bedfile = bedfile.sort_values("Start")

#export the file as a new bedfile
bedfile.to_csv(output, sep='\t', index=False, header=None)
