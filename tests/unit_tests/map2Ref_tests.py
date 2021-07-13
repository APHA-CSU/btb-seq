import unittest
from btb_tests import BtbTests

class Map2RefTests(BtbTests):
    ref_path = './references/Mycbovis-2122-97_LT708304.fas'

    def test_map2ref(self):
        """
            Asserts map2ref.bash completes on tinyreads without errors and produces a named bam file
        """        

        # Copy test data
        reads = self.copy_tinyreads(unzip=True)

        # Output Filenames
        output = self.temp_dirname + 'mapped.bam'

        # Pass case
        self.assertBashScript(0, ['./bin/map2Ref.bash', self.ref_path, reads[0], reads[1], output])
        self.assertFileExists(output)

        # Covert to SAM for visual checking
        self.bam_to_sam(output, self.temp_dirname+'mapped.sam')

if __name__ == '__main__':
    unittest.main()
