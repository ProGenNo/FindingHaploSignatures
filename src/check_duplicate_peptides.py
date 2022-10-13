import argparse
import pandas as pd

parser = argparse.ArgumentParser(
	description='Removes duplicate peptide sequences.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input TSV file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

peptides_df = pd.read_csv(args.input_file, sep='\t', header=0)

def group_rows(df):
    df['GeneID'] = ';'.join(df['GeneID'].tolist())
    df['ProteinID'] = ';'.join(df['ProteinID'].tolist())
    df['Position'] = df['Position'].apply(str)
    df['Position'] = ';'.join(df['Position'].tolist())
    df['type'] = ';'.join(df['type'].tolist())
    df['SNPs'] = ';'.join(df['SNPs'].tolist())
    df['matches_contaminant'] = any(df['matches_contaminant'].tolist())

    # aggregate the maximum number of SNPs for each gene (take a maximum for each gene ID)
    df['max_linked_SNPs'] = max([ max([ int(y) for y in x.split(',') if y.isdigit() ]) for x in df['possible_linked_SNPs'].tolist() ])
    df['possible_linked_SNPs'] = ';'.join([ str(max([ int(y) for y in x.split(',') ])) for x in df['possible_linked_SNPs'].tolist() ])

    # aggregate the peptide type
    type_aggregated = ""
    if 'canonical' in df['type_aggregated'].tolist():
        type_aggregated = "canonical"
    elif 'single_variant' in df['type_aggregated'].tolist():
        type_aggregated = "single_variant"
    elif 'multi_variant' in df['type_aggregated'].tolist():
        type_aggregated = 'multi_variant'

    df['type_aggregated'] = type_aggregated

    # store the peptide category
    category = ""
    if ';' in df['GeneID']:
        category = 'non_specific'
    elif ',' in df['ProteinID']:
        category = 'protein_specific'
    else:
        category = 'proteoform_specific'
    
    df['category'] = category

    return df

g = peptides_df.groupby(['Sequence']).apply(group_rows)
peptides_df['GeneID'] = g['GeneID']
peptides_df['ProteinID'] = g['ProteinID']
peptides_df['Position'] = g['Position']
peptides_df['type'] = g['type']
peptides_df['type_aggregated'] = g['type_aggregated']
peptides_df['category'] = g['category']
peptides_df['SNPs'] = g['SNPs']
peptides_df['matches_contaminant'] = g['matches_contaminant']
peptides_df['possible_linked_SNPs'] = g['possible_linked_SNPs']
peptides_df['max_linked_SNPs'] = g['max_linked_SNPs']
peptides_df.drop_duplicates(subset="Sequence", keep="first", inplace=True)

peptides_df.to_csv(args.output_file, sep='\t', header=True, index=False)
