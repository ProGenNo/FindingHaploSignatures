import pandas as pd
from common import get_protein_name_dict
import argparse

parser = argparse.ArgumentParser(
	description='Reads a PSM report, adds a column distinguishing non-specific (multigene) peptides.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input TSV file")

parser.add_argument("-f", dest="fasta_file", required=True,
                    help="input FASTA file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

print ('Reading', args.input_file)
pep_df = pd.read_csv(args.input_file, sep='\t', header=0)

print ('Reading', args.fasta_file)
name_dict, all_proteins = get_protein_name_dict(args.fasta_file)

def get_genes(proteins):
	protein_list = proteins.split(';')

	is_contaminant = not all([ (prot.startswith('ENSP') or prot.startswith('enshap')) for prot in protein_list ])
	if is_contaminant:
		return '-'

	ref_list = [ name_dict[prot].split(':')[0].split('.')[0] for prot in protein_list if prot.startswith('enshap') ]
	ref_list.extend([ prot.split('.', 1)[0] for prot in protein_list if prot.startswith('ENSP') ])
	ref_list = list(dict.fromkeys(ref_list))

	gene_list = []

	for prot in ref_list:
		gene_id = all_proteins[prot]['description'].split('gene:', 1)[1].split('.')[0]

		if gene_id not in gene_list:
			gene_list.append(gene_id)

	return ';'.join(gene_list)

pep_df['matching_genes'] = pep_df['Proteins'].apply(get_genes)

pep_df[['PSMId', 'Sequence', 'matching_genes']].to_csv(args.output_file, sep='\t', header=True, index=False)
