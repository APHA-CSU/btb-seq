set -e

### Ele's branches
# BaseQual - uses default base quality filter [13]
# MapQual - adds mapping quality filter as a parameter [30] 
# Ploidy - indicates haploid genome in bcftools call command [--ploidy 1] 
# VarQual - adds variant quality filter as a parameter and keeps snp-variants only [30] # also need to test params.VAR_QUAL = 20 
# ReadDepth - adds read depth filter (for high-quality reads) as a parameter [5] # also need to test params.MIN_READ_DEPTH = 7, 10
# ForRevAltRead - requires at least one read in each direction (F and R) in support of the alternate allele 
# AltProportion - adds filter that requires a minimum proportion of reads that support the alternate allele [0.75] # also need to test params.MIN_ALT_PROPORTION = 0.9, 0.95
# SNPwindow - removes SNPs that are within 10bp from each other [10] 
# RepeatMask - Masks sites that are within the repeat.bed file [references/Mycbovis-2122-97_LT708304.fas.rpt.regions.bed"] # also need to test 'references/Mbovis_mask_regions_adaptedPC2018.bed'

# Download each branch
root=/home/aaronfishman/ebs/pipeline-results/ofat-1/

python validator.py $root/branch/ --branch 

#
