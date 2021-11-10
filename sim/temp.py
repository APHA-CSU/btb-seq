import validator
import compare_snps
import glob
import os
import json
import pandas as pd

def analyse(root_path):
    """ Analyse results from an ofat run

        Parameters:
            root_path (str): Path to parent directory where 
                results from ofat trials are stored
    """

    # Determine path of each trial
    paths = glob.glob(root_path + '/*')

    # Initialise output data columns
    stats = []
    branch_names = []
    tps = []
    fns = []
    fps = []

    # Collect data from each trial
    for path in paths:
        branch_names.append(os.path.basename(path))

        pipeline_directory = glob.glob(path + '/btb-seq-results/Results*')[0] + '/'
        pipeline_snps = pipeline_directory + 'snpTables/simulated.tab'
        simulated_snps = path + '/simulated-genome/simulated.refseq2simseq.map.txt'
        
        try:        
            stats = compare_snps.analyse(simulated_snps, pipeline_snps)
            tps.append(stats['tp'])
            fns.append(stats['fn'])
            fps.append(stats['fp'])

        except Exception as e:
            print(e)

            tps.append('FAIL')
            fns.append('FAIL')
            fps.append('FAIL')

    

    return pd.DataFrame(data={
        "branch": branch_names,
        "fn": fns,
        "fp": fps,
        "tp": tps
    })

if __name__ == '__main__':
    analyse('/home/aaronfishman/temp-results/ofat-18/')
    quit()

    # Analyse Results
    # HACK: this could easily break if additioanl files are present
    pipeline_directory = glob.glob(btb_seq_results_path + '*')[0] + '/'
    pipeline_snps = pipeline_directory + 'snpTables/simulated.tab'
    stats = compare_snps.analyse(simulated_snps, pipeline_snps)


    pipeline_directory = glob.glob(btb_seq_results_path + '*')[0] + '/'


    pipeline_snps = '/home/aaronfishman/temp-results/ofat-18/AltProportion/btb-seq-results/Results_simulated-reads_09Nov21/snpTables/simulated.tab'
    simulated_snps = '/home/aaronfishman/temp-results/ofat-18/AltProportion/simulated-genome/simulated.refseq2simseq.map.txt'

    print(compare_snps.analyse(simulated_snps, pipeline_snps))
