#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import sys
import os

args=sys.argv

if len(args)<2:
    positions = "/home/prizamsandhu/mnt/fsx-036/Prizam/vcf_filtering_test/pos"
    mask = "/home/prizamsandhu/mnt/fsx-036/Prizam/vcf_filtering_test/DataDrivenMerge20.bed"
else:
    mask = sys.argv[1]
    positions =sys.argv[2]
    output =sys.argv[3]

f = open(positions).read().split('\n')

f = f[:-1]
pos = []

for x in f:
    pos.append(int(x))

#test if the samples has too many SNPs, if it does it is likey contaminated and will take too long to go through the 10 base SNP filter so output the old mask, will likey be caught by the other filters
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

bedfile = pd.DataFrame(columns=["Sample","Start","End"])
for x in chain:
    bedfile.loc[-1] = ["LT708304-Mycobacteriumbovis-AF2122-97",(x[0]-1),x[1]]
    bedfile.index = bedfile.index + 1
    bedfile = bedfile.sort_index()
    
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
