import argparse
import pandas as pd
import re
from multiprocessing import Pool
from tqdm import tqdm

parser = argparse.ArgumentParser(
	description='Removes duplicate peptide sequences.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input TSV file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

parser.add_argument("-t", dest="threads", type=int, required=False, default=5,
                    help="# threads to use")

parser.add_argument("-var", dest="variant_output_file", required=True,
                    help="output file")

args = parser.parse_args()

peptides_df = pd.read_csv(args.input_file, sep='\t', header=0)
all_lengths = peptides_df['Length'].drop_duplicates().tolist()

def group_rows(df):
    all_genes = ';'.join(df['GeneID'].tolist())
    df['GeneID'] = all_genes
    all_proteins = ';'.join(df['ProteinID'].tolist())
    df['ProteinID'] = all_proteins
    df['Position'] = df['Position'].apply(str)
    df['Position'] = ';'.join(df['Position'].tolist())
    all_types = ';'.join(df['type'].tolist())
    df['type'] = all_types
    df['SNPs'] = ';'.join(df['SNPs'].tolist())
    matches_contam = any(df['matches_contaminant'].tolist())
    df['matches_contaminant'] = matches_contam

    # aggregate the maximum number of SNPs for each gene (take a maximum for each gene ID)
    df['max_linked_SNPs'] = max([ max([ int(y) for y in x.split(',') if y.isdigit() ]) for x in df['possible_linked_SNPs'].tolist() ])
    df['possible_linked_SNPs'] = ';'.join([ str(max([ int(y) for y in x.split(',') ])) for x in df['possible_linked_SNPs'].tolist() ])

    # aggregate the peptide type
    all_types = re.split(r"[,;]", all_types)
    type_aggregated = ""
    if matches_contam:
        type_aggregated = "contaminant"
    elif 'canonical' in all_types:
        type_aggregated = "canonical"
    elif 'single_variant' in all_types:
        type_aggregated = "single_variant"
    elif 'multi_variant' in all_types:
        type_aggregated = 'multi_variant'

    df['type_aggregated'] = type_aggregated

    # store the peptide category
    category = ""
    if ';' in all_proteins:
        category = 'non_specific'
    elif ',' in all_proteins:
        category = 'protein_specific'
    else:
        category = 'proteoform_specific'
    
    df['category'] = category

    return df

def group_peptides(l):
    grouped_df = peptides_df[peptides_df['Length'] == l].groupby(['Sequence']).apply(group_rows)
    grouped_df.drop_duplicates(subset=['Sequence'], inplace=True, keep='first')
    return grouped_df

with Pool(min(args.threads, len(all_lengths))) as p:
    grouped_dfs = list(tqdm(p.imap_unordered(group_peptides, range(min(all_lengths), max(all_lengths) + 1)), total=(len(all_lengths))))
    print ("Writing output to", args.output_file)
    grouped_df = pd.concat(grouped_dfs)
    grouped_df.to_csv(args.output_file, sep='\t', header=True, index=False)

    print ("Writing variant peptides only to", args.variant_output_file)
    grouped_df[grouped_df['type_aggregated'].str.contains('variant')].to_csv(args.variant_output_file, sep='\t', header=True, index=False)
