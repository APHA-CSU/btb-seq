#!/usr/local/bin/python

""" unit_tests.py: 
        This file runs unit tests for the pipeline. 
        Each unit test should be quick to run (ie. milliseconds typically, but certainly less than second) 
        Tests that have a long run time should simplified, mocked, or be a .circleci job.
"""

import unittest

from deduplicate_tests import DeduplicateTests
from trim_tests import TrimTests
from mask_tests import MaskTests
from read_stats_tests import ReadStatsTests
from map2Ref_tests import Map2RefTests
from varCall_tests import VarCallTests

def main():
    # Test List
    tests = [
        DeduplicateTests,
        TrimTests,
        MaskTests,
        ReadStatsTests,
        Map2RefTests,
        VarCallTests
    ]

    # Load 'em up
    for test in tests:
        unittest.TestLoader().loadTestsFromTestCase(test)

    unittest.main()

if __name__ == '__main__':
    main()


