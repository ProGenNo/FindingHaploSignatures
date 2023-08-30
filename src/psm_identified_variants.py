import argparse
import pandas as pd


parser = argparse.ArgumentParser(
	description='Reads the PSM report file, creates a new file with an overview of identified amino acid variants.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input tab-separated file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

psm_df = pd.read_csv(args.input_file, sep='\t', header=0)

result_data = []
result_columns = ['ProteinVariant', 'Peptide', 'PeptideType', 'PSMId', 'Samples']

# store each found variant in each peptide as a separate entry, aggregate later
for index, row in psm_df.iterrows():
    coveredSNPs = row['CoveredSNPs'].split(';')
    peptide_type_global = row['PeptideType']
    for entry in coveredSNPs:
        if entry is '-':
            continue

        protein_stable_id = entry.split(':')[0]
        SNPs = entry.split(':')[1].split(',')
        peptide_type = ''

        if len(SNPs) > 1:
            peptide_type = 'multi-variant'
        elif len(SNPs) == 1:
            peptide_type = 'single-variant'
        else:
            print ('PSMID:', row['PSMId'], "no SNPs found in protein", protein_stable_id, entry)

        if (peptide_type != peptide_type_global):
            print ('PSMID:', row['PSMId'], "peptide downgraded -> SNPs not considered")
            continue

        for SNP in SNPs:
            result_data.append([protein_stable_id + ':' + SNP, row['Peptide'], peptide_type, row['PSMId'], row['sample_ID']])

result_df = pd.DataFrame(data=result_data, columns=result_columns)

def group_rows(df):
    df['Peptide'] = ';'.join(df['Peptide'].tolist())
    df['PeptideType'] = ';'.join(df['PeptideType'].tolist())
    df['PSMId'] = ';'.join(df['PSMId'].tolist())
    df['Samples'] = ';'.join(df['Samples'].tolist())
    return df

g = result_df.groupby(['ProteinVariant']).apply(group_rows)
result_df['Peptide'] = g['Peptide']
result_df['PeptideType'] = g['PeptideType']
result_df['PSMId'] = g['PSMId']
result_df['Samples'] = g['Samples']
result_df.drop_duplicates(subset="ProteinVariant", keep="first", inplace=True)

result_df.to_csv(args.output_file, sep='\t', header=True, index=False)
