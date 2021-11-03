import validator
import os

#
results_root = '/home/aaronfishman/ebs/pipeline-results/ofat-4/'
results_root = os.path.join(results_root, '')
repo_path = '/home/aaronfishman/temp2/btb-seq'
reference_path = '/home/aaronfishman/repos/btb-seq/sim/Mycobacterium_bovis_AF212297_LT78304.fa'

# First download them
branches = [
    "BaseQual",
    "MapQual",
    "Ploidy",
    "VarQual",
    "ReadDepth",
    "ForRevAltRead",
    "AltProportion",
    "SNPwindow",
    "RepeatMask"
]

def main():

    os.makedirs(results_root, exist_ok=False)

    for branch in branches:
        results_path = results_root + branch + '/'

        try:
            validator.performance_test(
                results_path, 
                repo_path, 
                reference_path, 
                branch=branch
            )
        except:
            print(f"***FAILED BRANCH: {branch}****", branch)

if __name__ == '__main__':
    main()