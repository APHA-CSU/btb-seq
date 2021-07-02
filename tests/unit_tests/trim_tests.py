import unittest
from btb_tests import BtbTests

class TrimTests(BtbTests):
    adapter_path = './references/adapter.fasta'

    def test_trim(self):
        """
            This introductory unit test asserts trim.bash completes on tinyreads without errors.
            And produces two fastq files
        """        

        # Copy test data
        reads = self.copy_tinyreads(unzip=True)

        # Output Filenames
        outputs = [
            self.temp_dirname + 'out1.fastq',
            self.temp_dirname + 'out2.fastq'
        ]

        # Pass case
        self.assertBashScript(0, ['./bin/trim.bash', self.adapter_path, reads[0], reads[1], outputs[0], outputs[1]])
        self.assertFileExists(outputs[0])
        self.assertFileExists(outputs[1])

        # Failure Case: adapter sequence file not found
        self.assertBashScript(1, ['./bin/trim.bash', './__does__/not/exist/', reads[0], reads[1], outputs[0], outputs[1]])

    def test_trimmomatric_trims(self):
        """
            Asserts that low quality reads are removed using a high quality and low quality examples
        """
        high_quality_sequence = "A"*40
        low_quality_sequence = "T"*40

        # Create Test Sequences
        high_quality = self.seq_record(high_quality_sequence, quality=93)
        low_quality = self.seq_record(low_quality_sequence, quality=0)

        read_1 = self.temp_dirname+"read_1.fastq"
        read_2 = self.temp_dirname+"read_2.fastq"

        output_1 = self.temp_dirname+"output_1.fastq"
        output_2 = self.temp_dirname+"output_2.fastq"

        self.write_fastq(read_1, [low_quality, high_quality])
        self.write_fastq(read_2, [low_quality, high_quality])
         
        # Run Trimmomatic
        self.assertBashScript(0, ['./bin/trim.bash', self.adapter_path, read_1, read_2, output_1, output_2])

        # Assert Output
        trimmed = self.read_fastq(output_1)
        self.assertEqual(trimmed, high_quality_sequence)

        trimmed = self.read_fastq(output_2)
        self.assertEqual(trimmed, high_quality_sequence)

if __name__ == '__main__':
    unittest.main()