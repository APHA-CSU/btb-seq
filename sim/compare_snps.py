import pandas as pd

from Bio import SeqIO

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

    simulated_pos -= masked_pos
    pipeline_pos -= masked_pos

    # TP - true positive -(the variant is in the simulated genome and correctly called by the pipeline)
    tp = len(simulated_pos.intersection(pipeline_pos))

    # FP (the pipeline calls a variant that is not in the simulated genome),
    fp = len(pipeline_pos - simulated_pos)

    # FN SNP calls (the variant is in the simulated genome but the pipeline does not call it).
    fn =  len(simulated_pos - pipeline_pos)

    # Compute Performance Stats
    # precision (positive predictive value) of each pipeline as TP/(TP + FP), 
    precision = tp / (tp + fp)

    # recall (sensitivity) as TP/(TP + FN)
    sensitivity = tp / (tp + fn)

    # miss rate as FN/(TP + FN)
    miss_rate = fn / (tp + fn)

    #  total number of errors (FP + FN) per million sequenced bases
    total_errors = fp + fn

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "sensitivity": sensitivity,
        "miss_rate": miss_rate,
        "total_errors": total_errors
    }

#TODO: This may not be required if we can get away with using bcftools/vcftools
#      for comparisons. Leaving this here for convenience in case those tools aren't suitable 
def load_consensus(path):
    """ Load a consensus file. Returns the first record in a fasta file a string """

    for seq_record in SeqIO.parse(path, "fasta"):
        return str(seq_record)
