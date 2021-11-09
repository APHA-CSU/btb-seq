import pandas as pd
import subprocess
from Bio import SeqIO
from io import StringIO

"""
Calculate performance stats from simulated data
"""

def masked_positions(mask_filepath='../references/Mycbovis-2122-97_LT708304.fas.rpt.regions'):
    mask = pd.read_csv(mask_filepath,
        delimiter='\t',
        skiprows=[0,1],
        header=None,
        names=["CHROM", "START", "END", "RPT"]
    )

    # TODO: This is off by one compared to the Emergency Port validation doc
    #    why is that?
    masked_pos = []
    for i, row in mask.iterrows():
        masked_pos.extend(list(range(row['START'], row['END']+1)))

    return masked_pos

def analyse(simulated_snps, pipeline_snps):
    """ Compare simulated SNPs data from simuG against btb-seq's snpTable.tab
        If adjust == True: applies the mask to simulated SNPs and pipeline SNPs
        Returns a dictionary of performance stats
    """

    # Load
    simulated = pd.read_csv(simulated_snps, delimiter='\t')
    pipeline = pd.read_csv(pipeline_snps, delimiter='\t')

    # Extract SNP positions
    simulated_pos = set(simulated['ref_start'].values)
    pipeline_pos = set(pipeline['POS'].values)
    masked_pos = set(masked_positions())

    simulated_pos_adjusted = simulated_pos - masked_pos
    pipeline_pos_adjusted = pipeline_pos - masked_pos

    # TP - true positive -(the variant is in the simulated genome and correctly called by the pipeline)
    tp = simulated_pos.intersection(pipeline_pos_adjusted)
    tp_rate = len(tp)

    # FP (the pipeline calls a variant that is not in the simulated genome),
    fp = pipeline_pos_adjusted - simulated_pos_adjusted
    fp_rate = len(fp)

    # FN SNP calls (the variant is in the simulated genome but the pipeline does not call it).
    fn = simulated_pos_adjusted - pipeline_pos_adjusted
    fn_rate =  len(fn)

    # TPs excluded 
    masked_tp = masked_pos.intersection(simulated_pos.intersection(pipeline_pos))
    masked_tp_rate = len(masked_tp)

    # FPs excluded
    masked_fp = masked_pos.intersection(simulated_pos - pipeline_pos)
    masked_fp_rate = len(masked_fp)

    # FNs excluded
    masked_fn = masked_pos.intersection(simulated_pos - pipeline_pos)
    masked_fn_rate = len(masked_fn)

    # Compute Performance Stats
    # precision (positive predictive value) of each pipeline as TP/(TP + FP), 
    precision = tp_rate / (tp_rate + fp_rate)

    # recall (sensitivity) as TP/(TP + FN)
    sensitivity = tp_rate / (tp_rate + fn_rate)

    # miss rate as FN/(TP + FN)
    miss_rate = fn_rate / (tp_rate + fn_rate)

    #  total number of errors (FP + FN) per million sequenced bases
    total_errors = fp_rate + fn_rate

    make_details_file('/home/nickpestell/pipeline-results/btb-seq-9/simulated-genome','/home/nickpestell/pipeline-results/btb-seq-9/btb-seq-results/Results_simulated-reads_09Nov21', fn)

    return {
        "TP": tp_rate,
        "FP": fp_rate,
        "FN": fn_rate,
        "masked TPs": masked_tp_rate,
        "masked FPs": masked_fp_rate,
        "masked FNs": masked_fn_rate,
        "precision": precision,
        "sensitivity": sensitivity,
        "miss_rate": miss_rate,
        "total_errors": total_errors
    }

def make_details_file(sim_path, pl_path, sites):
    simulated_genome = load_consensus(sim_path+'/simulated.simseq.genome.fa')
    pipeline_genome = load_consensus(pl_path+'/consensus/simulated.fas')
    pipeline_vcf = load_vcf(pl_path+'/vcf/simulated.vcf.gz')
    mask = set(masked_positions())
    with open(pl_path+'/details.txt','w') as details_file:
        for i in sites:              
            details_file.write("##POSITION: {}".format(i)+'\n')
            details_file.write('##SIMULATED GENOME:\n')
            details_file.write(simulated_genome[i-51:i-1].lower()+simulated_genome[i-1]+simulated_genome[i:i+50].lower()+'\n')
            details_file.write('##PIPELINE GENOME:\n')
            details_file.write(pipeline_genome[i-49:i+1].lower()+pipeline_genome[i+1]+pipeline_genome[i+2:i+52].lower()+'\n')
            details_file.write('##MASK:\n')
            mask_local = ''
            sites_local = set(range(i-50,i+51))
            for j in range(i-50,i+51):
                if j in mask.intersection(sites_local):
                    mask_local += 'M'
                else:
                    mask_local += ' '  
            details_file.write(mask_local+'\n')
            details_file.write("#PIPELINE VCF:\n")
            details_file.write(pipeline_vcf.loc[pipeline_vcf['POS'] == i].to_csv(sep='\t')+'\n')
            details_file.write("\n\n")


    x = 0

#TODO: This may not be required if we can get away with using bcftools/vcftools
#      for comparisons. Leaving this here for convenience in case those tools aren't suitable 
def load_consensus(path):
    """ Load a consensus file. Returns the first record in a fasta file a string """

    seq_record = SeqIO.read(path, "fasta")
    return str(seq_record.seq)

def load_vcf(path):
    returncode = subprocess.run(['gunzip',
                                 '-k',
                                 path]).returncode
    if returncode:
        raise Exception("""*****
            %s
            cmd failed with exit code %i
          *****""" % ('gunzip', returncode))
    with open(path[:-3], 'r') as vcf_file:
        vcf_data = [l for l in vcf_file if not l.startswith('##')]
    vcf_df = pd.read_csv(StringIO(''.join(vcf_data)), delimiter = '\t')

    return vcf_df
    
