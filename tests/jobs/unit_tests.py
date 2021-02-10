#!/usr/local/bin/python

""" unit_tests.py: 
        This file runs unit tests for the pipeline. 
        Each unit test should be quick to run (ie. milliseconds typically, but certainly less than second) 
        Tests that have a long run time should simplified, mocked, or be a .circleci job.
"""

import unittest
import subprocess
import tempfile
import os
import shutil
import glob

class TestPipeline(unittest.TestCase):
    def setUp(self):
        """
            Create new temporary directory
        """
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_dirname = self.temp_dir.__enter__() + '/'

        self.temp_dirname = './'

        print("Temporary Directory setup: " + self.temp_dirname)

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

    def test_deduplicate(self):
        """
            This introductory unit test asserts the deduplicate process completes on tinyreads without errors.
            It currently fails on line 23 because the depend_path arg is not supplied.
            In a future PR, I would like to fix this by putting FastUniq on the Linux search path.
        """

        pair_id = self.temp_dirname + 'tinyreads'

        for file in glob.glob(r'./tests/data/tinyreads/*'):
            shutil.copy(file, self.temp_dirname)

        return_code = subprocess.run(['bash', '-e', './bin/deduplicate.bash', pair_id]).returncode

        self.assertEqual(return_code, 0)

    def test_mask(self):
        pair_id = self.temp_dirname+'test'
        sam_filepath = './tests/data/tinymatch.sam'
        rpt_mask = './references/Mycbovis-2122-97_LT708304.fas.rpt.regions'

        # Convert to SAM to BAM
        with open(pair_id + '.mapped.sorted.bam', 'w') as f:
            subprocess.call(['samtools', 'view', '-S', '-b', sam_filepath], stdout=f)
       
        # Test it
        self.assertBashScript(0, ['./bin/mask.bash', pair_id, rpt_mask])

if __name__ == '__main__':
    unittest.main()
