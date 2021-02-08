#!/usr/local/bin/python

# Merge the csv output files and determine whether samples are M.bovis positive or not, based on predetemined 
# thresholds for number and proportion of M.bovis reads.

import pandas as pd
import argparse
from datetime import datetime, date

date_out = date.today().strftime('%d%b%y')

def combine(assigned_csv, bovis_csv, seq_run, read_threshold, abundance_threshold):
    #Read Assigned Clade csv and replace blank cells  with 'NA'
    assigned_df = pd.read_csv(assigned_csv)

    # Read BovPos csv and change ID to either 'Negative' or 'Mycobacterium bovis',
    # depending on number and relative abundance of M.bovis reads.
    # Thresholds set by comparison with output from spoligotyping.
    # Fill empty value cells with zero 
    qbovis_df = pd.read_csv(bovis_csv)
    qbovis_df.loc[(qbovis_df['TotalReads'] >= read_threshold) & (qbovis_df['Abundance'] >= abundance_threshold), 'ID'] = 'Mycobacterium bovis'
    qbovis_df.loc[(qbovis_df['TotalReads'] < read_threshold) & (qbovis_df['Abundance'] < abundance_threshold), 'ID'] = 'Negative'
    qbovis_df.loc[(qbovis_df['TotalReads'] < read_threshold) & (qbovis_df['Abundance'] > abundance_threshold), 'ID'] = 'Inconclusive'
    qbovis_df.loc[(qbovis_df['TotalReads'] > read_threshold) & (qbovis_df['Abundance'] < abundance_threshold), 'ID'] = 'Inconclusive'
    qbovis_df['TotalReads'].fillna('0', inplace = True)
    qbovis_df['Abundance'].fillna('0', inplace = True)
    qbovis_df['ID'].fillna('Negative', inplace = True)

    #Merge dataframes and fill ID with M bovis if Pass, then any remaining blank cells with 'NA'
    finalout_df = pd.merge(assigned_df, qbovis_df, on = 'Sample', how = 'outer')
    finalout_df['ID'].fillna('Mycobacterium bovis', inplace = True)
    finalout_df.fillna('NA', inplace = True)
    finalout_df.set_index('Sample', inplace = True)
    finalout_df.sort_index(inplace = True)

    #Write to csv
    finalout_df.to_csv(""+str(seq_run)+"_FinalOut_"+str(date_out)+".csv")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('assigned_csv', help='path to AssignedWGSClade.csv')
    parser.add_argument('bovis_csv', help='path to Bovis.csv')
    parser.add_argument('seq_run', help='Unique sequencer run number')
    parser.add_argument('read_threshold', nargs='?', type=int, default=500, help='threshold for number of M.bovis reads')
    parser.add_argument('abundance_threshold', nargs='?', type=int, default=1, help='threshold for M.bovis abundance')

    args = parser.parse_args()

    combine(**vars(args))
