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

    def tearDown(self):
        """
            Cleanup temporary directory.
        """
        self.temp_dir.__exit__(None, None, None)

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
        in_filepath = './tests/data/dummy.sam'

        # Convert to SAM
        cmd = ['samtools', 'view', '-S', '-b', in_filepath]
        
        print(cmd)

        with open(self.temp_dirname+'test.bam', 'w') as f:
            subprocess.call(cmd, stdout=f)

        # Test it

        pass

if __name__ == '__main__':
    unittest.main()
