from pipeline_test import TestPipeline
import subprocess

class MaskTests(TestPipeline):
    def test_mask(self):
        """
            Asserts mask.bash completes without errors when 
            the supplied sam file contains no regions of zero coverage
        """
        pair_id = self.temp_dirname + 'test'

        sam_filepath = './tests/data/tinymatch.sam'
        rpt_mask = './references/Mycbovis-2122-97_LT708304.fas.rpt.regions'

        # Convert to SAM to BAM
        with open(pair_id + '.mapped.sorted.bam', 'w') as f:
            subprocess.call(['samtools', 'view', '-S', '-b', sam_filepath], stdout=f)
       
        # Test
        self.assertBashScript(0, ['./bin/mask.bash', pair_id, rpt_mask])