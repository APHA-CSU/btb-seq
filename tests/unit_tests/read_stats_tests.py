import subprocess
import shutil
import glob

from btb_tests import BtbTests

class ReadStatsTests(BtbTests):
    def assertReadStats(self, reads_path, name):
        """
            Asserts read stats returns a 0 exit code for the supplied test reads
            TODO: Neaten test code when readStats.bash has input args for all input data
        """
        fastq_1 = reads_path+'_S1_R1_X.fastq.gz'
        fastq_2 = reads_path+'_S1_R2_X.fastq.gz'
        pair_id = self.temp_dirname + name

        # Copy over 
        shutil.copy(fastq_1, self.temp_dirname)
        shutil.copy(fastq_2, self.temp_dirname)

        # Unzip
        fastq_files = glob.glob(self.temp_dirname + '*.gz')
        subprocess.run(['gunzip', '-k'] + fastq_files, cwd=self.temp_dirname)

        fastq_file = fastq_files[0][:-3]

        subprocess.run(['ln', '-s', fastq_file, pair_id+'_uniq_R1.fastq'])
        subprocess.run(['ln', '-s', fastq_file, pair_id+'_uniq_R2.fastq'])
        shutil.copy(fastq_file, pair_id+'_trim_R1.fastq')

        self.sam_to_bam('./tests/data/tinymatch.sam', pair_id+'.mapped.sorted.bam')

        # Test
        self.assertBashScript(0, ['./bin/readStats.bash', pair_id])            

    def test_read_stats_tinyreads(self):
        self.assertReadStats('./tests/data/tinyreads/tinyreads', 'tinyreads')

    def test_read_stats_tinysra(self):
        self.assertReadStats('./tests/data/tinysra/tinysra', 'tinysra')
