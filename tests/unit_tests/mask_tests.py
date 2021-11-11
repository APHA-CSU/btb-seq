from btb_tests import BtbTests
import shutil
import unittest

class MaskTests(BtbTests):
    rpt_mask = './references/Mycbovis-2122-97_LT708304.fas.rpt.regions'

    def test_mask(self):
        """
            Asserts mask.bash completes without errors when 
            the supplied sam file contains no regions of zero coverage
        """
        # Copy data
        sam_filepath = self.temp_dirname+'tinysam.sam'
        bam_filepath = self.temp_dirname+'tinybam.bam'
        mask_filepath = self.temp_dirname+'masked.bed'
        shutil.copy('./tests/data/tinymatch.sam', sam_filepath)        
        
        # Convert to bam
        self.sam_to_bam(sam_filepath, bam_filepath)

        # Test
        self.assertBashScript(0, ['./bin/mask.bash', self.rpt_mask, bam_filepath, mask_filepath])
        self.assertFileExists(mask_filepath)

if __name__ == '__main__':
    unittest.main()
