from btb_tests import BtbTests
import shutil
import unittest

class Vcf2ConsensusTests(BtbTests):
    ref_path = './references/Mycbovis-2122-97_LT708304.fas'
    
    def test_vcf2consensus(self):
        """
            Asserts vcf2Consensus.bash completes without errors on a minimal example
        """
        min_read_depth = str(5)
        min_allele_freq = str(0.8)
        snp_gap = str(5)
        
        # Copy data
        bed_filepath = self.temp_dirname + 'tinybed.bed'
        vcf_filepath = self.temp_dirname + 'tinyvariants.vcf.gz'

        shutil.copy('./tests/data/tinybed.bed', bed_filepath) 
        shutil.copy('./tests/data/tinyvariants.vcf.gz', vcf_filepath) 

        # Outputs
        consensus_filepath = self.temp_dirname + 'consensus.fas'
        snps_filepath = self.temp_dirname + 'snps.tab'

        # Test
        self.assertBashScript(0, ['./bin/vcf2Consensus.bash', 
            self.ref_path, 
            bed_filepath,
            vcf_filepath,
            consensus_filepath,
            snps_filepath,
            min_read_depth,
            min_allele_freq,
            snp_gap,
            "test"
        ])

        self.assertFileExists(consensus_filepath)
        self.assertFileExists(snps_filepath)

if __name__ == '__main__':
    unittest.main()
