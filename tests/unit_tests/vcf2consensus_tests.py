from btb_tests import BtbTests
import shutil
import unittest

class Vcf2ConsensusTests(BtbTests):

    def test_vcf2consensus(self):
        """
            Asserts vcf2Consensus.bash completes without errors on a minimal example
        """
        # Copy data
        ref_filepath = self.temp_dirname + 'ref.fas'
        masked_filepath = self.temp_dirname + 'masked.bed'
        vcf_filepath = self.temp_dirname + 'edge-cases.vcf.gz'
        regions_filepath = self.temp_dirname + 'regions.bed'

        shutil.copy('./tests/data/tinyref.fas', ref_filepath) 
        shutil.copy('./tests/data/edge-cases.bed', masked_filepath) 
        shutil.copy('./tests/data/edge-cases.vcf.gz', vcf_filepath) 
        shutil.copy('./tests/data/edge-cases_regions.bed', regions_filepath) 

        self.index(vcf_filepath)

        # Outputs
        consensus_filepath = self.temp_dirname + 'consensus.fas'
        snps_filepath = self.temp_dirname + 'snps.tab'

        # Test
        self.assertBashScript(0, ['./bin/vcf2Consensus.bash', 
            ref_filepath, 
            masked_filepath,
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
