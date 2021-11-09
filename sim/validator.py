import os
import subprocess
import glob
import json
import sys
import argparse
import shutil

from compare_snps import analyse

DEFAULT_REFERENCE_PATH = './Mycobacterium_bovis_AF212297_LT78304.fa'

def run(cmd, *args, **kwargs):
    """ Run a command and assert that the process exits with a non-zero exit code.
        See python's subprocess.run command for args/kwargs

        Parameters:
            cmd (list): List of strings defining the command, see (subprocess.run in python docs)
            cwd (str): Set surr

        Returns:
            None
    """
    # TODO: store stdout to a file
    returncode = subprocess.run(cmd, *args, **kwargs).returncode

    if returncode:
        raise Exception("""*****
            %s
            cmd failed with exit code %i
          *****""" % (cmd, returncode))

def simulate_genome(reference_path, output_path, num_snps=16000):
    run([
        "simuG.pl",
        "-refseq", reference_path,
        "-snp_count", str(num_snps),
        "-prefix", output_path + "simulated"
    ])


def simulate_reads(
    genome_fasta,
    output_path,
    num_read_pairs=150000,
    read_length=150,
    seed=1,
    outer_distance=330,
    random_dna_probability=0.01,
    rate_of_mutations=0,
    indel_mutation_fraction=0,
    indel_extension_probability=0,
    per_base_error_rate="0" # TODO: default to Ele's reccomendation? 0.001-0.01
):
    output_prefix = output_path + "simulated"
    
    # How dwgsim chooses to name it's output fastq files
    dwgsim_read_1 = output_prefix + ".bwa.read1.fastq.gz"
    dwgsim_read_2 = output_prefix + ".bwa.read2.fastq.gz"

    # Output fastq filenames compatible with btb-seq
    read_1 = output_prefix + "_S1_R1_X.fastq.gz"
    read_2 = output_prefix + "_S1_R2_X.fastq.gz"

    run([
        "dwgsim",
        "-e", str(per_base_error_rate),
        "-E", str(per_base_error_rate),
        "-i",
        "-d", str(outer_distance),
        "-N", str(num_read_pairs),
        "-1", str(read_length),
        "-2", str(read_length),
        "-r", str(rate_of_mutations),
        "-R", str(indel_mutation_fraction),
        "-X", str(indel_extension_probability),
        "-y", str(random_dna_probability),
        "-H",
        "-z", str(seed),
        genome_fasta,
        output_prefix
    ])

    # Rename output fastq files
    os.rename(dwgsim_read_1, read_1)
    os.rename(dwgsim_read_2, read_2)


def btb_seq(btb_seq_directory, reads_directory, results_directory):
    run(["bash", "./btb-seq", reads_directory,
         results_directory], cwd=btb_seq_directory)


def performance_test(results_path, btb_seq_path, reference_path, exist_ok=False, branch=None):
    """ Runs a performance test against the pipeline

        Parameters:
            btb_seq_path (str): Path to btb-seq code is stored
            results_path (str): Output path to performance test results
            reference_path (str): Path to reference fasta
            exist_ok (bool): Whether or not to throw an error if a results directory already exists
            branch (str): Checkout git branch on the repo (default None)

        Returns:
            None
    """
    # Add trailing slash
    btb_seq_path = os.path.join(btb_seq_path, '')
    results_path = os.path.join(results_path, '')

    # Validate Input
    if os.path.isdir(results_path):
        raise Exception("Output results path already exists")
    os.makedirs(results_path)

    if not os.path.isdir(btb_seq_path):
        raise Exception("Pipeline code repository not found")

    # Output Directories
    simulated_genome_path = results_path + 'simulated-genome/'
    simulated_reads_path = results_path + 'simulated-reads/'
    btb_seq_backup_path = results_path + 'btb-seq/'
    btb_seq_results_path = results_path + 'btb-seq-results/'

    # Paths to simulated reference genome and simulated SNPs file
    fasta_path = simulated_genome_path + 'simulated.simseq.genome.fa'
    simulated_snps = simulated_genome_path + "simulated.refseq2simseq.map.txt"

    # TODO: handle dwgsim vcf files. Make sure we are taking into account variants it might generate

    # Create Output Directories
    os.makedirs(simulated_genome_path, exist_ok=exist_ok)
    os.makedirs(simulated_reads_path, exist_ok=exist_ok)
    os.makedirs(btb_seq_results_path, exist_ok=exist_ok)    

    # Backup btb-seq code
    # TODO: exclude the work/ subdirectory from this operation.
    #   This could potentially copy large amounts of data
    #   from the work/ directory nextflow generates
    shutil.copytree(btb_seq_path, btb_seq_backup_path)

    # TODO: Use nextflow's method of choosing github branches
    #    the method handles automatic pulling of github branches.
    if branch:
        checkout(btb_seq_backup_path, branch)

    # Run Simulation
    simulate_genome(reference_path, simulated_genome_path)
    simulate_reads(fasta_path, simulated_reads_path)
    btb_seq(btb_seq_backup_path, simulated_reads_path, btb_seq_results_path)

    # Analyse Results
    # HACK: this could easily break if additioanl files are present
    pipeline_directory = glob.glob(btb_seq_results_path + 'Results_simulated-reads_*')[0] + '/'
    pipeline_snps = pipeline_directory + 'snpTables/simulated.tab'
    stats = analyse(simulated_snps, pipeline_snps)

    # Write output
    with open(results_path + "stats.json", "w") as file:
        file.write(json.dumps(stats, indent=4))

def checkout(repo_path, branch):
    run(["git", "checkout", str(branch)], cwd=repo_path)

def main():
    # Parse
    parser = argparse.ArgumentParser(
        description="Performance test btb-seq code")
    parser.add_argument("results", help="path to performance test results")
    parser.add_argument("--btb_seq", default="../", help="path to btb-seq code")
    parser.add_argument("--ref", "-r", help="path to reference fasta", default=DEFAULT_REFERENCE_PATH)
    parser.add_argument("--branch", help="path to reference fasta", default=None)

    args = parser.parse_args(sys.argv[1:])

    # Run
    performance_test(args.results, args.btb_seq, args.ref, args.branch)

if __name__ == '__main__':
    main()
