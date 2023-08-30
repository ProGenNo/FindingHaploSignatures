import argparse
import pandas as pd
import matplotlib.pyplot as plt

def sort_types(zipped):
    if zipped[0] == 'decoy':
        return 999999999
    return -zipped[1]

parser = argparse.ArgumentParser(
	description='Reads the PSM report file, creates violin plot for confidence measures by peptide type.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input tab-separated file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

parser.add_argument("-val_field", dest="value_field", required=True,
                    help="value field")

parser.add_argument("-thr_field", dest="threshold_field", required=True,
                    help="threshold field")

parser.add_argument("-thr_value", dest="threshold_value", required=True,
                    help="threshold value", type=float)
                    
parser.add_argument("-title", dest="figure_title", required=True,
                    help="figure title")

parser.add_argument("-ylabel", dest="y_label", required=True,
                    help="y axis label")

args = parser.parse_args()

print ("Reading", args.input_file)
psm_df = pd.read_csv(args.input_file, sep='\t', header=0)
psm_df = psm_df[(psm_df[args.threshold_field] < args.threshold_value) & (~psm_df['PeptideType'].str.match('contaminant')) & (~psm_df['PeptideType'].str.match('has_stop'))]

uniq_types = list(pd.unique(psm_df['PeptideType'].values.ravel()))
#uniq_types = ['canonical', 'protein-coding-noncanonical', 'contaminant']

# create the array to be plotted
fig_dataset = []
n_samples = []
for peptype in uniq_types:
    subset = psm_df[psm_df['PeptideType'] == peptype]
    fig_dataset.append(subset[args.value_field].values)
    n_samples.append(len(subset))

zipped = list(zip(uniq_types, n_samples, fig_dataset))
zipped.sort(key=sort_types)
uniq_types, n_samples, fig_dataset = zip(*zipped)

fig, axes = plt.subplots()

violin_parts = axes.violinplot(dataset = fig_dataset)

decoy_idx = uniq_types.index('decoy')

violin_parts['bodies'][decoy_idx].set_facecolor('red')
violin_parts['cbars'].set_color('grey')
violin_parts['cmins'].set_color('grey')
violin_parts['cmaxes'].set_color('grey')

axes.set_title(args.figure_title)

x_ticks = [ peptype + '\n#PSMs = {:,.0f}'.format(n_samples[i]) for i,peptype in enumerate(uniq_types) ]

axes.set_ylabel(args.y_label)
axes.set_xticks([i+1 for i in range(len(uniq_types))])
axes.set_xticklabels(x_ticks, rotation=45, fontsize=10)
fig.tight_layout()

plt.savefig(args.output_file)
#plt.show()

