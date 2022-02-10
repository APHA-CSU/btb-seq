from btb_tests import BtbTests
import shutil
import unittest

class MaskTests(BtbTests):
    allsites = './references/All-sites.bed'

    def test_mask(self):
        """
            Asserts mask.bash completes without errors when 
            the supplied sam file contains no regions of zero coverage
        """
        # Copy data
        rpt_mask_filepath = self.temp_dirname+'rpt_mask.bed'
        vcf_filepath = self.temp_dirname+'edge-cases.vcf.gz'
        mask_filepath = self.temp_dirname+'masked.bed'
        regions_filepath = self.temp_dirname+'regions.bed'
        shutil.copy('./tests/data/edge-cases.vcf.gz', vcf_filepath)        
        shutil.copy('./tests/data/tinymask.bed', rpt_mask_filepath)        
        
        # Test
        self.assertBashScript(0, ['./bin/mask.bash', 
            rpt_mask_filepath, 
            vcf_filepath, 
            mask_filepath, 
            regions_filepath, 
            self.allsites, 
            str(5), 
            str(0.8)])
        self.assertFileExists(regions_filepath)

if __name__ == '__main__':
    unittest.main()
