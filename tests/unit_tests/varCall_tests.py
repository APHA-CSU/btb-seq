import shutil
import unittest
from btb_tests import BtbTests

class VarCallTests(BtbTests):
    ref_path = './references/Mycbovis-2122-97_LT708304.fas'

    def test_varCall(self):
        """
            Asserts varCall.bash completes on tinymatch.sam without errors and produces a vcf file
        """
        sam_path = self.temp_dirname + '/tinymatch.sam'
        bam_path = self.temp_dirname + '/tinymatch.bam'
        output = self.temp_dirname + 'variants.vcf.gz'

        # Copy test data
        shutil.copy2('./tests/data/tinymatch.sam', sam_path)

        # Convert to BAM
        self.sam_to_bam(sam_path, bam_path)

        # Pass case
        self.assertBashScript(0, ['./bin/varCall.bash', self.ref_path, bam_path, output])
        self.assertFileExists(output)

        # Unzip
        self.unzip(output)

if __name__ == '__main__':
    unittest.main()
