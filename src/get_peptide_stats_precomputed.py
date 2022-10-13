from common import get_protein_name_dict
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
	description='Reads a precomputed peptide stats file, gives summary statistics on protein coverage by canonical/non-canonical peptides.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input TSV file, peptide stats by protein")

parser.add_argument("-f", dest="fasta_file", required=True,
                    help="input FASTA file")

args = parser.parse_args()

name_dict, all_proteins = get_protein_name_dict(args.fasta_file)

pep_df = pd.read_csv(args.input_file, sep='\t', header=0)

total_aa = [0, 0, 0, 0]
total_var = [0, 0, 0]

for index, row in pep_df.iterrows():
    coverage_str = row['coverage']

    for region_str in coverage_str.split(';'):
        region = list(map(lambda x: int(x), region_str.split('_')))
        region_len = region[1] - region[0]
        region_type = region[2] + 1
        total_aa[region_type] += region_len

    total_var[0] += int(row['SNP stats'].split('single_variant:')[1].split(',')[0])
    total_var[1] += int(row['SNP stats'].split('multi_variant:')[1].split(',')[0])
    total_var[2] += int(row['SNP stats'].split('both:')[1].split(',')[0])


total_aa_sum = 0

for protein in all_proteins.values():
    if (protein['accession'].startswith('ENSP')):
        total_aa_sum += len(protein['sequence'])

total_var_sum = sum(total_var)

print ("Proteome length:", total_aa_sum)
print ("Canonical proteome: %d AAs - %.2f %%,\npossible single-variant peptides: %d AAs - %.2f %%,\npossible multi_variant peptides: %d AAs - %.2f %%,\nsequences not matching to peptides: %d AAs - %.2f %%" % (total_aa[1], (total_aa[1] / total_aa_sum) * 100, total_aa[2], (total_aa[2] / total_aa_sum) * 100, total_aa[3], (total_aa[3] / total_aa_sum) * 100, total_aa[0], (total_aa[0] / total_aa_sum) * 100))
print ()
print ("Total number of variants:", total_var_sum)
print ("Variants covered in single-var. pept.: %d - %.2f %%; \nvariants covered in both single-var. and multi-var. pept.: %d - %.2f %%; \nvariants covered in multi_variant peptides: %d - %.2f %%" % (total_var[0], (total_var[0] / total_var_sum) * 100, total_var[1], (total_var[1] / total_var_sum) * 100, total_var[2], (total_var[2] / total_var_sum) * 100))
