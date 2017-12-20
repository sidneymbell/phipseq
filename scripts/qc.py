import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import argparse
from scipy.stats import linregress
from collections import defaultdict

######### Input ########

## Read in data; use the ID to index individual oligos
args = argparse.ArgumentParser()
args.add_argument('-input', type=str, default='../data/2017-09-28/annotatedCounts.xlsx',
                    help='Path to input counts excel sheet')
args.add_argument('-out_path', type=str, default='../data/2017-09-28/')
args.add_argument('-samples', type=str, default='../data/2017-09-28/samples.tsv')
args = args.parse_args()

counts = pd.read_excel(args.input, index_col=0)
out_path = args.out_path

#########  Initial cleanup & normalization #########
 # We know these animals were previously vaccinated with HIV antigen; drop those oligos
hiv_oligos = counts.loc[counts['Virus_Strain'].str.contains('HIV')]
counts.drop(hiv_oligos.index.values, axis=0, inplace=True)

# For convenience, separate out the metadata, input counts, background counts and metadata.
def proportions(df):
    ''' For each column, divide each element by the sum of the column (column sums to 1)'''
    xsum=df.sum(0)
    df = df.div(xsum, axis='columns')
    return df

input_cols = [ c for c in counts.columns.values if 'input' in c.lower() ]
beads_cols = [ c for c in counts.columns.values if 'beads' in c.lower() ]
metadata_cols = ['Virus_Strain', 'Start_to_End_nt', 'Peptide_sequence']
sample_cols = [c for c in counts if not any([c in input_cols, c in beads_cols, c in metadata_cols])]

metadata = counts[metadata_cols]
# Standardize each column to sum to 1
values = proportions(counts[[c for c in counts.columns.values if c not in metadata_cols]])
input_ctrls = values[input_cols]
beads_ctrls = values[beads_cols]

######  Find technical and biological replicates, plot direct comparisons #########
sample_sera_map = { i.split()[0]: i.split()[1] for i in open(args.samples, 'r').readlines() }
biological_replicates = defaultdict(list)
for sample, sera in sample_sera_map.items():
    all_sample_reps = [c for c in sample_cols if sample.upper() in c.upper()]
    biological_replicates[sera] += all_sample_reps
biological_replicates = dict(biological_replicates)

technical_replicates = defaultdict(list)
for serum in sample_cols:
    technical_replicates[serum.rsplit('-', 1)[0]].append(serum) # Find replicates
technical_replicates = dict(technical_replicates) # Turn off defaultdict behavior

sns.set(style='whitegrid', font_scale = 1.3, palette='pastel') ## Make all of our plots prettier
def compare_replicates(df, columns, title, fname):
    ''' Plot sanity checks for technical replicates '''
    if columns:
        if len(columns) == 1:
            return
        replicates = df[columns]
        replicates.fillna(0, inplace=True)
    else:
        replicates = df

    def plot_comparison(x,y, **kwargs):
        scatter = plt.plot(x,y, 'o', alpha=0.4)
        try:
            r_2 = linregress(x,y)[2]
            scatter[0].axes.text(0,0, 'R^2 = %.2f'%(r_2))
        except:
            pass
        return scatter

    g = sns.PairGrid(replicates, diag_sharey=False)
    g.map_diag(sns.violinplot)
    g.map_offdiag(plot_comparison, )

    g.fig.suptitle(title, va='bottom')
    plt.tight_layout()
    plt.savefig(out_path+'/figs/'+fname, bbox_inches='tight')
    plt.close()

for serum, tech_reps in technical_replicates.items():
    compare_replicates(values, tech_reps, serum, serum+'_tech_reps.png')

for serum, bio_reps in biological_replicates.items():
    compare_replicates(values, bio_reps, serum, serum+'_bio_reps.png')

compare_replicates(values, input_cols + beads_cols, 'Input + Beads', 'ctrl_reps.png')
#######  Tidy up the metadata a bit  ########
# Tidy start and end coordinates --> integers
metadata['start'], metadata['end'] = metadata['Start_to_End_nt'].str.split('to', 1).str
metadata['start'] = metadata['start'].map(lambda x: int(x))
metadata['end'] = metadata['end'].map(lambda x: int(x.split('.')[0])) ## TODO
metadata.drop('Start_to_End_nt', inplace=True, axis=1)

# Tidy up virus and strain names
def parse_strains(virusstrain):
    # e.g., 'DENV3_BR-BID-V2403-2008.DENV3_Mozambique1985'
    # The ONNV sequences overlap with the CHIKV sequences; for now, we'll omit it; this should be revisited. TODO

    names = [s for s in virusstrain.split('.') if 'ONNV' not in s] # ['DENV3_BR-BID-V2403-2008', 'DENV3_Mozambique1985']
    virus = [s.split('_', 1)[0] for s in names] # ['DENV3']

    if len(set(virus)) != 1:
        virus, strains = np.nan, np.nan

    else:
        virus = virus[0] # 'DENV3'
        strains = [s.split(virus+'_', 1)[1].replace('-', '').replace('_', '').upper() for s in names if s != '']
        # ['BRBIDV24032008', 'MOZAMBIQUE1985']

    return pd.Series({'virus': virus, 'strains':strains})

new_names = [parse_strains(v) for v in metadata['Virus_Strain']]
metadata['virus'] = [n['virus'] for n in new_names]
metadata['strains'] = [n['strains'] for n in new_names]
metadata.rename(columns={'Peptide_sequence':'sequence'}, inplace=True)
metadata.drop('Virus_Strain', axis=1, inplace=True)
metadata.dropna(how='any', inplace=True, axis=(0, 1))

#####   Write to file   #####
values = values.join(metadata, how='inner')
values.to_csv(out_path+'proportions.csv')

print 'Put plots to compare technical and biological replicates to %s.\nAdd the name of any sample that needs to be dropped (one per line) to: %sdrop.txt'%(out_path+'/figs/', out_path)
