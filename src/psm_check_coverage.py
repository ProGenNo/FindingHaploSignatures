import argparse
import re
import pandas as pd
from common import read_fasta

parser = argparse.ArgumentParser(
	description='Reads the PSM report file, creates a new file with information about covered redions for each candidate protein.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input tab-separated file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

parser.add_argument("-f", dest="fasta_file", required=True,
                    help="fasta file")                    

args = parser.parse_args()

all_proteins = read_fasta(args.fasta_file)

psm_df = pd.read_csv(args.input_file, sep='\t', header=0)
psm_df = psm_df[(psm_df['PeptideType'] != 'decoy') & (psm_df['PeptideType'] != 'contaminant') & (psm_df['PeptideType'] != 'has_stop')]

summary_data = []
summary_columns = ['ProteinID', 'ProteinHaplotype', 'CoveredRegions']

coverage = {}

for index, row in psm_df.iterrows():
    protein_accessions = re.split(r',|;', row['Proteins'])   # all the accessions of candidate proteins
    peptide_length = len(row['Sequence'])                 # length of the peptide sequence
    peptide_starts = list(map(lambda x: int(x), row['Position'].replace(',', ';').split(';')))    # positions of the peptide within respective candidate protein

    for prot_idx, acc in enumerate(protein_accessions):
        protein_id = '-'
        protein_haplotype = acc
        if ('enshap' in acc and acc in all_proteins):
            protein_id = all_proteins[acc]['description'].split('ref:')[1].split('.')[0]
            protein_haplotype = all_proteins[acc]['description'].split('haplotype:')[1].split()[0]
        elif (acc.startswith('ENSP')):
            protein_id = acc.split('.')[0]
            protein_haplotype = 'REF'

        region = str(peptide_starts[prot_idx]) + '-' + str(peptide_starts[prot_idx] + peptide_length)

        if acc in coverage:
            if region in coverage[acc]['regions']:
                coverage[acc]['regions'][region] += 1
            else:
                coverage[acc]['regions'][region] = 1
        else:
            coverage[acc] = { 'protein_id': protein_id, 'regions': { region: 1 }, 'haplotype' : protein_haplotype }

for acc in coverage:
    regions = []
    for region in coverage[acc]['regions']:
        regions.append( region + "(" + str(coverage[acc]['regions'][region]) + ")")
    regions.sort(key=lambda x: int(x.split('-')[0]))
    summary_data.append([coverage[acc]['protein_id'], coverage[acc]['haplotype'], ",".join(regions)])

summary_df = pd.DataFrame(data=summary_data, columns=summary_columns)
summary_df.sort_values(by=['ProteinID'], inplace=True)
summary_df.to_csv(args.output_file, sep='\t', header=True, index=False)
