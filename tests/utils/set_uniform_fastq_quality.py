#!/usr/local/bin/python

""" set_uniform_fastq_quality.py: 
        Copy a fastq to a new file with uniform base quality
"""

import sys
import argparse

from Bio import SeqIO

def set_uniform_fastq_quality(quality, filepath_in, filepath_out):
    """
        Copy a fastq file to a new file with uniform base quality
    """
    # Parse
    original_fastq = SeqIO.parse(filepath_in, "fastq")

    # Obliterate quality
    new_fastq = []
    for i, line in enumerate(original_fastq):
        num_reads = len(line._per_letter_annotations["phred_quality"])
        new_qualities = [quality]*num_reads
        line._per_letter_annotations["phred_quality"] = new_qualities
        
        new_fastq.append(line)

    # Write
    SeqIO.write(new_fastq, filepath_out, "fastq")

if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(description="Copy a fastq to a new file with the worst possible quality")
    parser.add_argument('quality', type=int, help="Phred quality score that all bases are set to")
    parser.add_argument('filepath_in', help="Filepath to original fastq file")
    parser.add_argument('filepath_out', help="Output filepath")

    args = parser.parse_args()

    # Run
    set_uniform_fastq_quality(**vars(args))
