import validator
import os
import sys
import argparse

"""
    Run the performance benchmarking tool against a number of git branches
"""

# Initial list of branches that we are testing
DEFAULT_BRANCHES = [
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

def ofat(btb_seq_path, results_path, branches=DEFAULT_BRANCHES):
    """ Runs a performance test against the pipeline

        Parameters:
            btb_seq_path (str): Path to btb-seq code is stored
            results_path (str): Output path to performance test results
            branches (list): List of strings for each git branch to test
    """
    # Add trailing slash
    btb_seq_path = os.path.join(btb_seq_path, '')
    results_path = os.path.join(results_path, '')     
    
    # Default reference for sims
    reference_path = btb_seq_path + './sim/Mycobacterium_bovis_AF212297_LT78304.fa'

    # Prepare output directory
    os.makedirs(results_path, exist_ok=False)

    # Benchmark the branches
    for branch in branches:
        branch_results_path = results_path + branch + '/'

        print('branch_results_path', branch_results_path)
        print('btb_seq_path', btb_seq_path)
        print('reference_path', reference_path)

        try:
            validator.performance_test(
                branch_results_path, 
                btb_seq_path, 
                reference_path, 
                branch=branch
            )
        except Exception as e:
            print(e)
            print(f"***FAILED BRANCH: {branch}****", branch)

if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(description="Run the performance benchmarking tool against a number of git branches")

    parser.add_argument("btb_seq", help="path to btb-seq code")
    parser.add_argument("results", help="path to results directory")

    args = parser.parse_args(sys.argv[1:])

    # Run
    ofat(args.btb_seq, args.results)