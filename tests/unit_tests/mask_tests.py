from btb_tests import BtbTests
import shutil
import unittest

class MaskTests(BtbTests):
    rpt_mask = './references/Mycbovis-2122-97_LT708304.fas.rpt.regions'
    allsites = './references/All-sites.bed'

    def test_mask(self):
        """
            Asserts mask.bash completes without errors when 
            the supplied sam file contains no regions of zero coverage
        """
        # Copy data
        vcf_filepath = self.temp_dirname+'tinyvariants.vcf.gz'
        mask_filepath = self.temp_dirname+'masked.bed'
        regions_filepath = self.temp_dirname+'regions.bed'
        shutil.copy('./tests/data/tinyvariants.vcf.gz', vcf_filepath)        
        
        # Test
        self.assertBashScript(0, ['./bin/mask.bash', 
            self.rpt_mask, 
            vcf_filepath, 
            mask_filepath, 
            regions_filepath, 
            self.allsites, 
            str(5), 
            str(0.8)])
        self.assertFileExists(mask_filepath)

if __name__ == '__main__':
    unittest.main()
