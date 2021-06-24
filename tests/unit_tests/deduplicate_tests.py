from btb_tests import BtbTests
import glob
import shutil

class DeduplicateTests(BtbTests):
    def test_deduplicate(self):
        """
            This introductory unit test asserts the deduplicate process completes on tinyreads without errors.
        """

        # Copy test data
        pair_id = self.temp_dirname + 'tinyreads'

        for file in glob.glob(r'./tests/data/tinyreads/*'):
            shutil.copy(file, self.temp_dirname)

        # Test the script
        self.assertBashScript(0, ['./bin/deduplicate.bash', pair_id])