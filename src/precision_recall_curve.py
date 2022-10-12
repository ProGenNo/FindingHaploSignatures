import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
        description='Plot the precision-recall curve based on q-values.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input tab-separated file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

parser.add_argument("-thr_field", dest="threshold_field", required=True,
                    help="threshold field")

parser.add_argument("-thr_max", dest="max_threshold", required=True,
                    help="maximum threshold value", type=float)

parser.add_argument("-xlabel", dest="x_label", required=True,
                    help="x axis label")

args = parser.parse_args()

print ("Reading", args.input_file)
psm_df = pd.read_csv(args.input_file, sep='\t', header=0)
psm_df = psm_df[(psm_df['PeptideType'] != 'decoy')]
target_psm_count = len(psm_df)
psm_df = psm_df[(psm_df[args.threshold_field] <= args.max_threshold)]

def get_count(group_df):
	group_df['count'] = len(group_df)
	return group_df

g = psm_df.groupby(args.threshold_field).apply(get_count)
g = g[[args.threshold_field, 'count']].sort_values(by=args.threshold_field).drop_duplicates()

g['cumsum'] = g['count'].cumsum()
g['percent_accepted'] = g['cumsum'] / target_psm_count * 100

plt.figure()
plt.plot(g[args.threshold_field].tolist(),g['percent_accepted'].tolist())
plt.xlabel(args.x_label)
plt.ylabel("% accepted target PSMs")
plt.grid()

plt.savefig(args.output_file)


