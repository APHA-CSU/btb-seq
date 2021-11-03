import glob
import os
import json
import pandas as pd

root = '/home/aaronfishman/ebs/pipeline-results/ofat-4/*'

# Loop over each subfolder
paths = glob.glob(root)

fns = []
fps = []
tps = []
branch_names = []

for path in paths:
    branch = os.path.basename(path)

    branch_names.append(branch)

    stats_path = path+'/stats.json'

    if not os.path.exists(stats_path):
        fn = 'FAIL'
        fp = 'FAIL'
        tp = 'FAIL'

    else:
        # Load the data 
        with open(stats_path, 'r') as file:
            data = json.load(file)

        fn = data['fn']
        fp = data['fp']
        tp = data['tp']

    fns.append(fn)
    fps.append(fp)
    tps.append(tp)

df = pd.DataFrame(data={
    "branch": branch_names,
    "fn": fns,
    "fp": fps,
    "tp": tps
})

print(df)