from btb_tests import BtbTests
import shutil
import unittest

class Vcf2ConsensusTests(BtbTests):
    ref_path = './references/Mycbovis-2122-97_LT708304.fas'

    def test_vcf2consensus(self):
        """
            Asserts vcf2Consensus.bash completes without errors on a minimal example
        """
        # Copy data
        mask_filepath = self.temp_dirname + 'tinybed.bed'
        vcf_filepath = self.temp_dirname + 'tinyvariants.vcf.gz'
        regions_filepath = self.temp_dirname + 'regions.bed'

        shutil.copy('./tests/data/tinybed.bed', mask_filepath) 
        shutil.copy('./tests/data/tinyvariants.vcf.gz', vcf_filepath) 
        shutil.copy('./tests/data/tinybed.bed', regions_filepath) 

        self.index(vcf_filepath)

        # Outputs
        consensus_filepath = self.temp_dirname + 'consensus.fas'
        snps_filepath = self.temp_dirname + 'snps.tab'

        # Test
        self.assertBashScript(0, ['./bin/vcf2Consensus.bash', 
            self.ref_path, 
            mask_filepath,
            regions_filepath,
            vcf_filepath,
            consensus_filepath,
            snps_filepath,
            "test.bcf"
        ])
        self.assertFileExists(consensus_filepath)
        self.assertFileExists(snps_filepath)


if __name__ == '__main__':
    unittest.main()
