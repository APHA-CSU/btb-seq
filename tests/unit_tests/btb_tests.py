import unittest
import subprocess
import tempfile
import os
import glob
import shutil
import subprocess

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class BtbTests(unittest.TestCase):
    """
        Base class for btb-seq unit tests
        Contains tools that make writing unit tests for nextflow processes easier
    """
    def setUp(self):
        """
            Create new temporary directory
        """
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_dirname = self.temp_dir.__enter__() + '/'

        print('Temp Dir: ' + self.temp_dirname)

    def tearDown(self):
        """
            Cleanup temporary directory.
        """
        self.temp_dir.__exit__(None, None, None)

    def assertBashScript(self, expected_exit_code, cmd):
        """
            Runs a bash script with exit on error (set -e) and throws an exception if the exit code doesn't match the expected 
            
            Inputs: 
                error_code <int>
                cmd <list> The command as a list of strings. eg ['hello-world.bash', '-name', 'Aaron']
        """

        actual_exit_code = subprocess.run(['bash', '-e'] + cmd).returncode
        self.assertEqual(expected_exit_code, actual_exit_code)

    def assertFileExists(self, filepath):
        """
            Raises an exception if the file does not exist
        """
        if not os.path.exists(filepath):
            raise Exception(f"File does not exist: {filepath}")

    def sam_to_bam(self, sam_filepath, bam_filepath):
        # Convert to SAM to BAM
        with open(bam_filepath, 'w') as f:
            subprocess.call(['samtools', 'view', '-S', '-b', sam_filepath], stdout=f)
    
    def bam_to_sam(self, bam_filepath, sam_filepath):
        # Convert to BAM to SAM
        proc = subprocess.run(['samtools', 'view', '-h', '-o', sam_filepath, bam_filepath])
        self.assertEquals(0, proc.returncode)

    def copy_tinyreads(self, unzip=False):
        """
            Copies tinyreads to the temporary directory that tests run in
        """
        for read in glob.glob(r'./tests/data/tinyreads/*.fastq.gz'):
            shutil.copy2(read, self.temp_dirname)

        reads = glob.glob(self.temp_dirname + '*.fastq.gz')

        # Unzip
        if unzip:
            for read in reads:
                self.unzip(read)

            reads = [read[:-3] for read in reads]

        # Return path to files
        return reads

    def unzip(self, path):
        """ Unzip a .gz file inplace """
        proc = subprocess.run(['gunzip', path])
        self.assertFalse(proc.returncode)

    def write_fastq(self, filepath, seq_records):
        """
            Write a fastq file to the filepath

            filepath:  (str) filepath of output fastq file
            sequences: (list/str) a list of strings that represent each sequence. Can also provide just a string if there is only one sequence
            quality:   (int) uniform phred33 quality score for each base
        """
        with open(filepath, "w") as file:
            Bio.SeqIO.write(seq_records, file, "fastq")

    def read_fastq(self, filepath):
        """
            Reads a fastq file and returns the sequence as a string
            Throws an exception if there is more than one record
        """
        return str(SeqIO.read(filepath, "fastq").seq)

    def seq_record(self, seq_str, quality=93):
        """
            Makes a new Seq Record with uniform quality from a string of ATCG's 
        """
        seq = SeqRecord(Seq(seq_str))
        seq.letter_annotations["phred_quality"] = [quality]*len(seq)
        return seq
