import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm

parser = argparse.ArgumentParser(
        description='Make the stacked bar chart of relative peptide categories by variation type.')

parser.add_argument("-i", dest="unique_list", required=True,
                    help="input unique peptides TSV file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

pep_df = pd.read_csv(args.unique_list, sep='\t', header=0)

type_dfs = [ pep_df[pep_df['type_aggregated'] == 'canonical'], pep_df[pep_df['type_aggregated'] == 'single_variant'], pep_df[pep_df['type_aggregated'] == 'multi_variant'] ]

labels = ['non-specific', 'protein-specific', 'proteoform-specific']
types = ['canonical', 'single-vairant', 'multi-variant']
stack_data = []
bar_width = 0.75

print('quantities: non-specific, protein-specific, proteoform-specific')

for i,type_df in enumerate(type_dfs):
	cat1 = len(type_df[type_df['GeneID'].str.contains(';')])
	cat2 = len(type_df[(~type_df['GeneID'].str.contains(';')) & type_df['ProteinID'].str.contains(',')])
	cat3 = len(type_df[(~type_df['GeneID'].str.contains(';')) & (~type_df['ProteinID'].str.contains(','))])

	print(types[i], cat1, cat2, cat3)

	stack_data.append([cat1 / len(type_df) * 100, cat2 / len(type_df) * 100, cat3 / len(type_df) * 100])

fig, ax = plt.subplots()

nonsp_vals = [ stack_data[0][0], stack_data[1][0], stack_data[2][0] ]
prot_vals = [ stack_data[0][1], stack_data[1][1], stack_data[2][1] ]
proteof_vals = [ stack_data[0][2], stack_data[1][2], stack_data[2][2] ]
combined_vals = [ stack_data[0][1] + stack_data[0][2], stack_data[1][1] + stack_data[1][2], stack_data[2][1] + stack_data[2][2] ]

ax.bar(labels, proteof_vals, bar_width, color=plt.cm.inferno(0.25), label='proteoform-specific')
ax.bar(labels, prot_vals, bar_width, bottom = proteof_vals, color=plt.cm.inferno(0.5), label='protein-specific')
ax.bar(labels, nonsp_vals, bar_width, bottom = combined_vals, color=plt.cm.inferno(0.8), label='non-specific')

x_ticks = ['canonical\nn = {:,.0f}'.format(len(type_dfs[0])), 'single-variant\nn = {:,.0f}'.format(len(type_dfs[1])), 'multi-variant\nn = {:,.0f}'.format(len(type_dfs[2]))]
ax.set_xticklabels(x_ticks, rotation=45, fontsize=10)

ax.set_ylabel('Share of Peptides [%]')
ax.set_xlabel('Peptide Type')
#ax.set_title('')
ax.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=1)

fig.tight_layout()

plt.savefig(args.output_file, dpi=600)

'''
print('quantities: non-specific, protein-specific, proteoform-specific')
print('Canonical' + str(stack_data[0]))
print('Single-variant' + str(stack_data[1]))
print('Multi-variant' + str(stack_data[2]))
'''




