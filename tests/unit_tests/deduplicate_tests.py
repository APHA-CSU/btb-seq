import unittest
import glob
import shutil

from btb_tests import BtbTests

class DeduplicateTests(BtbTests):
    def test_deduplicate(self):
        """
            This introductory unit test asserts the deduplicate process completes on tinyreads without errors.
            And produces two fastq files
        """

        # Copy test data
        reads = glob.glob(r'./tests/data/tinyreads/*') 
        outputs = [self.temp_dirname+'1.txt', self.temp_dirname+'2.txt']

        # Test the script
        self.assertBashScript(0, ['./bin/deduplicate.bash', reads[0], reads[1], outputs[0], outputs[1]])
        self.assertFileExists(outputs[0])
        self.assertFileExists(outputs[1])

        # TODO: Assert that duplicates are actually removed


if __name__ == '__main__':
    unittest.main()