#!/usr/local/bin/python

""" obliterate_fastq_quality.py: 
        Copy a fastq to a new file with the worst possible quality
"""

import sys
import argparse

from Bio import SeqIO

def obliterate_quality(filepath_in, filepath_out):
    # Parse
    original_fastq = SeqIO.parse(filepath_in, "fastq")

    # Obliterate quality
    badq_fastq = []
    for line in original_fastq:
        new_qualities = [0]*len(line._per_letter_annotations["phred_quality"])
        line._per_letter_annotations["phred_quality"] = new_qualities
        badq_fastq.append(line)

    # Write
    SeqIO.write(badq_fastq, filepath_out, "fastq")

if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(description="Copy a fastq to a new file with the worst possible quality")
    parser.add_argument('filepath_in', help="Filepath to original fastq file")
    parser.add_argument('filepath_out', help="Output filepath")

    args = parser.parse_args()

    # Run
    obliterate_quality(**vars(args))
