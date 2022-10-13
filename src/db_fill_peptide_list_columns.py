import argparse
import pandas as pd
from multiprocessing import Pool
from tqdm import tqdm

parser = argparse.ArgumentParser(
	description='Fill missing columns in the peptide list aggregated by protein.')

parser.add_argument("-nonuniq", dest="list_by_protein", required=True,
                    help="input TSV file")

parser.add_argument("-uniq", dest="unique_list", required=True,
                    help="input TSV file")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="# threads to use")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

list_by_protein = pd.read_csv(args.list_by_protein, sep='\t', header=0)

unique_list = pd.read_csv(args.unique_list, sep='\t', header=0)
unique_list.set_index('Sequence', inplace=True)

def get_category_contaminant(seq):
    sequence_row = unique_list.loc[seq]
    return [ sequence_row['category'], sequence_row['matches_contaminant'] ]

def process_row(i):
    row = list_by_protein.iloc[i]
    result = get_category_contaminant(row['Sequence'])

    return result

with Pool(args.threads) as p:
    #perm = np.random.permutation(total_protein_count)
    all_results = list(tqdm(p.imap(process_row, range(0, len(list_by_protein))), total=len(list_by_protein)))

    categories = [ elem[0] for elem in all_results ]
    is_contam = [ elem[1] for elem in all_results ]

    list_by_protein['category'] = categories
    list_by_protein['matches_contaminant'] = is_contam

    list_by_protein.to_csv(args.output_file, sep='\t', header=True, index=False)