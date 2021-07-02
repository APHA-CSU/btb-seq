import unittest
from btb_tests import BtbTests

class TrimTests(BtbTests):
    def test_trim(self):
        """
            This introductory unit test asserts trim.bash completes on tinyreads without errors.
            And produces two fastq files
        """
        adapter_path = './references/adapter.fasta'

        # Copy test data
        reads = self.copy_tinyreads(unzip=True)

        # Output Filenames
        outputs = [
            self.temp_dirname + 'out1.fastq',
            self.temp_dirname + 'out2.fastq'
        ]

        # Pass case
        self.assertBashScript(0, ['./bin/trim.bash', adapter_path, reads[0], reads[1], outputs[0], outputs[1]])
        self.assertFileExists(outputs[0])
        self.assertFileExists(outputs[1])

        # Failure Case: adapter sequence file not found
        self.assertBashScript(1, ['./bin/trim.bash', './__does__/not/exist/', reads[0], reads[1], outputs[0], outputs[1]])

        # TODO: Assert that adapters are indeed trimmed


if __name__ == '__main__':
    unittest.main()