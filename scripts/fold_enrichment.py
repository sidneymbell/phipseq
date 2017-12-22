import pandas as pd
import numpy as np
import argparse
from collections import defaultdict
from scipy.stats import linregress


#####   Parse Input   #####
args = argparse.ArgumentParser()
args.add_argument('-input', type=str, default='../data/2017-09-28/proportions.csv',
                    help='Path to input oligo proportions')
args.add_argument('-out_path', type=str, default='../data/2017-09-28/')
args.add_argument('-samples', type=str, default='../data/2017-09-28/samples.tsv',
                    help='Tab-delimited file of sampleID    seraID, e.g.: NHP-1    ZIKV')
args.add_argument('-drop', type=str, default='../data/2017-09-28/drop.txt',
                    help='Text file with which samples to drop, one sample per line.')
args.add_argument('--aggregate_technical', type=bool, default=True,
                    help='Aggregate by technical replicates? Default True')
args.add_argument('--aggregate_biological', type=bool, default=False,
                    help='Aggregate by biological replicate? Default False')
args.add_argument('-control_col', type=str, default='beads',
                    help='Column name to use as baseline when calculating fold enrichment')
args = args.parse_args()

out_path = args.out_path
sample_sera_map = { i.split()[0]: i.split()[1] for i in
                    open(args.samples, 'r').readlines() }   # { 'NHP-3': 'ZIKV' }

drop_samples = [line.strip() for line in open(args.drop, 'r')]  #[ 'NHP-3-1', 'NHP-3-2']
control_col = args.control_col # 'input'

# pd.DataFrame(index=oligoID, columns=sampleID, values=proportions of reads in each column assigned to each oligo)
proportions = pd.read_csv(args.input, index_col=0)
# drop bad replicates or superfluous samples as specified
proportions.drop(drop_samples, inplace=True, axis=1, errors='ignore')

# all columns with 'input' in the name
input_cols = [c for c in proportions.columns.values if 'input' in c.lower()]
# all columns with 'beads' in the name
beads_cols = [c for c in proportions.columns.values if 'beads' in c.lower()]

# nonnumerical columns
metadata_cols = ['virus', 'start', 'end', 'Peptide_sequence', 'start', 'end', 'strains']
metadata = proportions[metadata_cols]

# all non-metadata columns
proportions = proportions[[c for c in proportions.columns.values if c not in metadata_cols]]

#####   Aggregate replicates    #####
def aggregate(df, reps, name):
    if len(reps) == 1:
        df[name] = df[reps[0]]
    else:
        df[name] = df[reps].mean(axis=1) # mean of each row --> new column of aggregated values
    df.drop(reps, inplace=True, axis=1, errors='ignore') # drop original columns (yes, the axis designator switches between these two methods which is incredibly annoying but is correct I promise)

if args.aggregate_technical:
    technical_replicates = defaultdict(list)
    # {'NHP-3-2ng': ['NHP-3-2ng-1', 'NHP-3-2ng-2']}

    for serum in proportions.columns.values: # Find replicates like ['NHP-3-1', 'NHP-3-2']
        if serum in beads_cols:
            name = 'beads'
        elif serum in input_cols:
            name = 'input'
        else:
            name = serum.rsplit('-', 1)[0]
        technical_replicates[name].append(serum)

    for serum, reps in technical_replicates.items():
        aggregate(proportions, reps, serum)
    # now data looks like pd.DataFrame(columns=['NHP-1', 'input', 'beads', 'NHP-2', ...])

if args.aggregate_biological:
    assert args.aggregate_technical, 'ERROR: cannot aggregate biological replicates without first aggregating technical replicates'
    biological_replicates = defaultdict(list) # {'ZIKV': ['NHP-6-2ng', 'NHP-6-10ng', 'NHP-3']}

    for sample, sera in sample_sera_map.items():
        if sample in proportions.columns.values:
            biological_replicates[sera].append(sample)

    for serum, reps in biological_replicates.items():
        aggregate(proportions, reps, serum)

#####   Convert all values to fold enrichment over control  #####
def fold_enrichment(row):
    ctrl = row[control_col]
    def enrichment_value(i):
        if ctrl == 0.:
            ctrl = proportions[control_col].mean() # for the rare 0 values in the control, plug in the mean control value
        return float(i) / float(ctrl)
    return row.map(enrichment_value)

mean_control_val = proportions[control_col].mean() # avoid dividing by 0
filled_control_col = proportions[control_col].replace(0., mean_control_val)
enrichment = proportions.divide(filled_control_col, axis=0)

annotated_enrichment = enrichment.join(metadata) # reattach the metadata
annotated_enrichment.to_csv(out_path+'enrichment.csv') # write to file
