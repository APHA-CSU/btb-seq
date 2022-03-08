from btb_tests import BtbTests
import shutil
import unittest

class MaskTests(BtbTests):
    allsites = './references/All-sites.bed'
    ref_masked_filepath = './tests/data/edge-cases.bed'

    def test_mask(self):
        """
            Asserts mask.bash completes without errors when 
            the supplied sam file contains no regions of zero coverage
        """
        # Copy data
        rpt_mask_filepath = self.temp_dirname+'rpt_mask.bed'
        vcf_filepath = self.temp_dirname+'edge-cases.vcf'
        vcf_gz_filepath = self.temp_dirname+'edge-cases.vcf.gz'
        masked_filepath = self.temp_dirname+'masked.bed'
        regions_filepath = self.temp_dirname+'regions.bed'
        shutil.copy('./tests/data/edge-cases.vcf', vcf_filepath)        
        shutil.copy('./tests/data/tinymask.bed', rpt_mask_filepath)        
        
        self.gzip(vcf_filepath, vcf_gz_filepath)

        # Test
        self.assertBashScript(0, ['./bin/mask.bash', 
            rpt_mask_filepath, 
            vcf_gz_filepath, 
            masked_filepath, 
            regions_filepath, 
            self.allsites, 
            str(8), 
            str(0.8),
            str(0.5)])
        self.assertFileExists(regions_filepath)
        with open(masked_filepath) as test_f, open(self.ref_masked_filepath) as ref_f:
            self.assertListEqual(list(test_f), list(ref_f))


if __name__ == '__main__':
    unittest.main()
