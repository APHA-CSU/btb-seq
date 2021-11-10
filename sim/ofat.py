import os
import sys
import argparse
import glob
import json

import pandas as pd
import validator


"""
    Run the performance benchmarking tool against a multiple git branches
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

def ofat(btb_seq_path, results_path, reference_path, branches=DEFAULT_BRANCHES):
    """ Runs a performance test against the pipeline

        Parameters:
            btb_seq_path (str): Path to btb-seq code is stored
            results_path (str): Output path to performance test results
            branches (list): List of strings for each git branch to test
    """
    # Add trailing slash
    btb_seq_path = os.path.join(btb_seq_path, '')
    results_path = os.path.join(results_path, '')     
    
    # TODO: set the reference file relative to the 
    #   This is awkward to do atm because many of the branches
    #   we are testing are not up to date with the latest sim code

    # Prepare output directory
    os.makedirs(results_path, exist_ok=False)

    # Benchmark the branches
    for branch in branches:
        branch_results_path = results_path + branch + '/'

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

    # Analyse results
    return analyse(results_path)

def analyse(root_path):
    """ Analyse results from an ofat run

        Parameters:
            root_path (str): Path to parent directory where 
                results from ofat trials are stored
    """

    # Determine path of each trial
    paths = glob.glob(root_path + '/*')

    # Initialise output data columns
    fns = []
    fps = []
    tps = []
    branch_names = []

    # Collect data from each trial
    for path in paths:
        # Initialised trial data
        branch = os.path.basename(path)
        stats_path = path + '/stats.json'

        fn = 'FAIL'
        fp = 'FAIL'
        tp = 'FAIL'

        # Log 
        if os.path.exists(stats_path):
            with open(stats_path, 'r') as file:
                data = json.load(file)

            fn = data['fn']
            fp = data['fp']
            tp = data['tp']

        # Add to output
        fns.append(fn)
        fps.append(fp)
        tps.append(tp)
        branch_names.append(branch)

    return pd.DataFrame(data={
        "branch": branch_names,
        "fn": fns,
        "fp": fps,
        "tp": tps
    })


if __name__ == '__main__':
    # Parse
    parser = argparse.ArgumentParser(description="Run the performance benchmarking tool against a number of git branches")

    parser.add_argument("btb_seq", help="path to btb-seq code")
    parser.add_argument("results", help="path to results directory")
    parser.add_argument("reference", help="path to reference file")

    args = parser.parse_args(sys.argv[1:])

    # Run
    df = ofat(args.btb_seq, args.results, args.reference)
    print("OFAT run completed, results:")
    print(df)