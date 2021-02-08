#!/usr/local/bin/python

# merge the csv output files and determine whether samples are M.bovis positive or not, based on predetemined 
# thresholds for number and proportion of M.bovis reads.

import pandas as pd
import sys, os
from datetime import datetime, date

In1 = sys.argv[1]
In2 = sys.argv[2]
SeqRun = sys.argv[3]

DateOut = date.today().strftime('%d%b%y')

#Read Assigned Cluster csv and replace blank cells  with 'NA'
Assigned = pd.read_csv(In1)

# Read BovPos csv and change ID to either 'Negative' or 'Mycobacterium bovis',
# depending on number and relative abundance of M.bovis reads.
# Thresholds set by comparison with output from spoligotyping.
# Fill empty value cells with zero 
Qbovis = pd.read_csv(In2)
Qbovis.loc[(Qbovis['TotalReads'] >= 500) & (Qbovis['Abundance'] >= 1), 'ID'] = 'Mycobacterium bovis'
Qbovis.loc[(Qbovis['TotalReads'] < 500) & (Qbovis['Abundance'] < 1), 'ID'] = 'Negative'
Qbovis.loc[(Qbovis['TotalReads'] < 500) & (Qbovis['Abundance'] > 1), 'ID'] = 'Inconclusive'
Qbovis.loc[(Qbovis['TotalReads'] > 500) & (Qbovis['Abundance'] < 1), 'ID'] = 'Inconclusive'
Qbovis['TotalReads'].fillna('0', inplace = True)
Qbovis['Abundance'].fillna('0', inplace = True)
Qbovis['ID'].fillna('Negative', inplace = True)

#Merge dataframes and fill ID with M bovis if Pass, then any remaining blank cells with 'NA'
FinalOut = pd.merge(Assigned, Qbovis, on = 'Sample', how = 'outer')
FinalOut['ID'].fillna('Mycobacterium bovis', inplace = True)
FinalOut.fillna('NA', inplace = True)
FinalOut.set_index('Sample', inplace = True)
FinalOut.sort_index(inplace = True)

#Write to csv
FinalOut.to_csv(""+str(SeqRun)+"_FinalOut_"+str(DateOut)+".csv")