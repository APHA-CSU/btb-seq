#!/usr/local/bin/python

""" unit_tests.py: 
		This file runs unit tests for the pipeline. 
		Each unit test should be quick to run (ie. less than second). 
		Tests that have a long run time should be an .circleci job.
"""

import unittest
import subprocess

class TestPipeline(unittest.TestCase):
    def test_deduplicate(self):
        """
            This introductory unit test asserts the deduplicate process completes on tinyreads without errors.
            It currently fails on line 23 because the depend_path arg is not supplied.
            In a future PR, I would like to fix this by putting FastUniq on the Linux search path.
        """
        pair_id = 'tests/data/tinyreads/tinyreads'

        return_code = subprocess.run(['bash', '-e', './bin/deduplicate.bash', pair_id]).returncode

        self.assertEqual(return_code, 0)

if __name__ == '__main__':
    unittest.main()
