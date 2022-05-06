#!/usr/bin/env python3

# Merge the csv output files and determine whether samples are M.bovis positive or not, based on predetemined 
# thresholds for number and proportion of M.bovis reads.

import pandas as pd
import argparse
import getpass
import git
import os
from datetime import datetime, date

def combine(assigned_csv, bovis_csv, seq_run, commitId, read_threshold, abundance_threshold):

    date_out = date.today().strftime('%d%b%y')
    user = getpass.getuser()
    scriptpath = os.path.dirname(os.path.abspath(__file__))
    repo = git.Repo(scriptpath, search_parent_directories=True)
    #commit = repo.head.object.__str__()

    #Read Assigned Clade csv and replace blank cells  with 'NA'
    assigned_df = pd.read_csv(assigned_csv)
    assignedround_df = assigned_df.round(2)

    # Read BovPos csv and change ID to either 'Negative' or 'Mycobacterium bovis',
    # depending on number and relative abundance of M.bovis reads.
    # Thresholds set by comparison with output from spoligotyping.
    # Fill empty value cells with zero 
    # TODO: Potential that not all scenarios are covered by the options below, so need to add default ID outcome
    bovis_df = pd.read_csv(bovis_csv)
    qbovis_df = bovis_df.round(2)
    qbovis_df.loc[(qbovis_df['TotalReads'] >= read_threshold) & (qbovis_df['Abundance'] >= abundance_threshold), 'ID'] = 'Mycobacterium bovis'
    qbovis_df.loc[(qbovis_df['TotalReads'] < read_threshold) & (qbovis_df['Abundance'] < abundance_threshold), 'ID'] = 'Negative'
    qbovis_df.loc[(qbovis_df['TotalReads'] < read_threshold) & (qbovis_df['Abundance'] > abundance_threshold), 'ID'] = 'Inconclusive'
    qbovis_df.loc[(qbovis_df['TotalReads'] > read_threshold) & (qbovis_df['Abundance'] < abundance_threshold), 'ID'] = 'Inconclusive'
    qbovis_df['TotalReads'].fillna('0', inplace = True)
    qbovis_df['Abundance'].fillna('0', inplace = True)
    qbovis_df['ID'].fillna('Negative', inplace = True)

    #Merge dataframes fill with appropriate Mycobacterium ID (Other, microti, bovis), then any remaining blank cells with 'NA'
    finalout_df = pd.merge(assignedround_df, qbovis_df, on = 'Sample', how = 'outer')
    finalout_df.loc[(finalout_df['group'] == 'nonbTB' ) | (finalout_df['group'] == 'MicPin' ) | (finalout_df['group'] == 'Pinnipedii' ), 'ID' ] = 'Other Mycobacteria'
    finalout_df.loc[(finalout_df['group'] == 'Microti' ), 'ID' ] = 'Mycobacterium microti'
    finalout_df['ID'].fillna('Mycobacterium bovis', inplace = True)
    finalout_df.fillna('NA', inplace = True)
    finalout_df.set_index('Sample', inplace = True)
    finalout_df.sort_index(inplace = True)

    #Write to csv
    finalout_df.to_csv("{}_FinalOut_{}.csv".format(seq_run, date_out))

    #Append log info
    with open("{}_FinalOut_{}.csv".format(seq_run, date_out), "a") as outFile:
        outFile.write("# Operator: " +user +"\n" +"# BovTB-nf commit: " +commitId)
        outFile.close

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('assigned_csv', help='path to AssignedWGSClade.csv')
    parser.add_argument('bovis_csv', help='path to Bovis.csv')
    parser.add_argument('seq_run', help='Unique sequencer run number')
    parser.add_argument('commitId', help='Nextflow capture of git commit')
    parser.add_argument('--read_threshold', type=int, default=500, help='threshold for number of M.bovis reads')
    parser.add_argument('--abundance_threshold', type=int, default=1, help='threshold for M.bovis abundance')

    args = parser.parse_args()

    combine(**vars(args))
