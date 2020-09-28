#!/usr/local/bin/python

""" set_uniform_fastq_quality.py: 
        Overwrite a fastq with uniform base quality
"""

import sys
import argparse

from Bio import SeqIO

def set_uniform_fastq_quality(quality, filepath):
    """
        Overwrite a fastq file to a new file with uniform base quality
    """
    # Parse
    original_fastq = SeqIO.parse(filepath, "fastq")

    # Obliterate quality
    new_fastq = []
    for i, line in enumerate(original_fastq):
        num_reads = len(line._per_letter_annotations["phred_quality"])
        new_qualities = [quality]*num_reads
        line._per_letter_annotations["phred_quality"] = new_qualities
        
        new_fastq.append(line)

    # Write
    SeqIO.write(new_fastq, filepath, "fastq")


if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(description="Overwrite a fastq to with a uniform quality")
    parser.add_argument('quality', type=int, help="Phred quality score that all bases are set to")
    parser.add_argument('filepath', help="Filepath to fastq file")

    args = parser.parse_args()

    # Run
    set_uniform_fastq_quality(**vars(args))
