import validator
import os
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
    repo_path = os.path.join(repo_path, '')
    results_path = os.path.join(results_path, '')     
    
    # Default reference for sims
    reference_path = repo_path + './sim/Mycobacterium_bovis_AF212297_LT78304.fa'

    # Prepare output directory
    os.makedirs(results_path, exist_ok=False)

    # Run against branches
    for branch in branches:
        branch_results_path = results_path + branch + '/'

        try:
            validator.performance_test(
                branch_results_path, 
                repo_path, 
                reference_path, 
                branch=branch
            )
        except:
            print(f"***FAILED BRANCH: {branch}****", branch)

if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(description="Run the performance benchmarking tool against a number of git branches")

    parser.add_argument("results", help="path to results directory")
    parser.add_argument("btb_seq", help="path to btb-seq code")

    args = parser.parse_args(sys.argv[1:])

    # Run
    ofat(args.results, args.btb_seq)