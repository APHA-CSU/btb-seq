from btb_tests import BtbTests
import shutil
import unittest

class Vcf2ConsensusTests(BtbTests):
    ref_consensus_filepath = './tests/data/edge-cases.fas'

    def test_vcf2consensus(self):
        """
            Asserts vcf2Consensus.bash completes without errors on a minimal example
        """
        # Copy data
        ref_filepath = self.temp_dirname + 'ref.fas'
        masked_filepath = self.temp_dirname + 'masked.bed'
        vcf_filepath = self.temp_dirname + 'edge-cases.vcf'
        vcf_gz_filepath = self.temp_dirname+'edge-cases.vcf.gz'
        regions_filepath = self.temp_dirname + 'regions.bed'
        bcf_filepath = self.temp_dirname + 'edge-cases.bcf'

        shutil.copy('./tests/data/tinyref.fas', ref_filepath) 
        shutil.copy('./tests/data/edge-cases.bed', masked_filepath) 
        shutil.copy('./tests/data/edge-cases.vcf', vcf_filepath) 
        shutil.copy('./tests/data/edge-cases_regions.bed', regions_filepath) 

        self.gzip(vcf_filepath, vcf_gz_filepath)
        self.index(vcf_gz_filepath)

        # Outputs
        consensus_filepath = self.temp_dirname + 'consensus.fas'
        unmasked_consensus_filepath = self.temp_dirname + 'unmasked_consensus.fas'
        snps_filepath = self.temp_dirname + 'snps.tab'

        # Test
        self.assertBashScript(0, ['./bin/vcf2Consensus.bash', 
            ref_filepath, 
            masked_filepath,
            regions_filepath,
            vcf_gz_filepath,
            consensus_filepath,
            unmasked_consensus_filepath,
            snps_filepath,
            bcf_filepath,
            str(0.8)
        ])
        self.assertFileExists(snps_filepath)
        with open(consensus_filepath) as test_f, open(self.ref_consensus_filepath) as ref_f:
            self.assertMultiLineEqual(test_f.read(), ref_f.read())


if __name__ == '__main__':
    unittest.main()
