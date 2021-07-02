import unittest
import subprocess
import tempfile
import os
import glob
import shutil
import subprocess

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
                proc = subprocess.run(['gunzip', read])
                self.assertFalse(proc.returncode)

            reads = [read[:-3] for read in reads]

        # Return path to files
        return reads
