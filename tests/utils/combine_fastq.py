#!/usr/local/bin/python

""" combine_fastq.py: 
        Creates a new fastq file by taking random subsets from two fastq files
"""

import argparse
import random
import gzip

def combine_fastq(
    total_num_reads, 
    prop_a_reads, 
    filepath_a, 
    filepath_b, 
    filepath_combined,
    seed
    ):
    random.seed(seed)

    num_rand_a_reads = int(total_num_reads*prop_a_reads)
    num_rand_b_reads = int(total_num_reads - num_rand_a_reads)

    # Setup output file
    with gzip.open(filepath_a, 'r') as file_a, \
        gzip.open(filepath_b, 'r') as file_b, \
        gzip.open(filepath_combined, 'w') as file_combined:

        subsample(num_rand_a_reads, file_a, file_combined)
        subsample(num_rand_b_reads, file_b, file_combined)

def subsample(num_rand_reads, file_in, file_out):
    """
        Randomly subsample fastq file_in, writing to file_out
    """
    # Get the number of records
    file_in.seek(0)
    num_reads=int(len(file_in.readlines())/4)
    file_in.seek(0)

    # Subsample
    record_nums = sorted(random.sample(range(num_reads), num_rand_reads))
    rec_num = -1
    for rr in record_nums:
        # Loop until next record
        while rec_num < rr:
            rec_num += 1
            for i in range(4): file_in.readline()

        # Write the random record
        for i in range(4):
            file_out.write(file_in.readline())
        rec_num += 1

if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(description="Overwrite a fastq to with a uniform quality")
    parser.add_argument('total_num_reads', type=float, help="Number of reads to take from first fastq file")
    parser.add_argument('prop_a_reads', type=float, help="Number of reads to take from second fastq file")
    parser.add_argument('filepath_a', help="First fastq filepath")
    parser.add_argument('filepath_b', help="Second fastq filepath")
    parser.add_argument('filepath_combined', help="Output combined fastq filepath")
    parser.add_argument('--seed', default=None, type=int, help="seed to initialise random number generator")

    args = parser.parse_args()

    # Run
    combine_fastq(**vars(args))
