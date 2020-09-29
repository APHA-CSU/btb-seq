#!/usr/local/bin/python

""" set_uniform_fastq_quality.py: 
        Overwrite a fastq with uniform base quality
"""

import argparse

from Bio import SeqIO

def set_uniform_fastq_quality(quality, filepath):
    """
        Overwrite a fastq file to a new file with uniform base quality
    """
    # Parse
    original_fastq = SeqIO.parse(filepath, "fastq")

    # Set quality
    new_fastq = []
    for i, line in enumerate(original_fastq):
        num_reads = len(line._per_letter_annotations["phred_quality"])
        line._per_letter_annotations["phred_quality"] = [quality]*num_reads
        
        new_fastq.append(line)

    # Write
    SeqIO.write(new_fastq, filepath, "fastq")


if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(description="Overwrite a fastq to with a uniform quality")
    parser.add_argument('quality', type=int, help="Phred quality score that all bases are set to")
    parser.add_argument('filepath', nargs="+", help="Filepath(s) to fastq file")

    args = parser.parse_args()

    # Run
    for filepath in args.filepath:
        set_uniform_fastq_quality(args.quality, filepath)
