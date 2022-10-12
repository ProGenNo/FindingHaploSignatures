import argparse
import pandas as pd

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

parser.add_argument("-thr_field", dest="threshold_field", required=True,
                    help="threshold field")

parser.add_argument("-thr_value", dest="threshold_value", required=True,
                    help="threshold value", type=float)

args = parser.parse_args()

psm_df = pd.read_csv(args.input_file, sep='\t', header=0)
result_df = psm_df[(psm_df[args.threshold_field] < args.threshold_value) & (~psm_df['PeptideType'].str.match('contaminant')) & (~psm_df['PeptideType'].str.match('has_stop'))]

result_df.to_csv(args.output_file, sep='\t', header=True, index=False)