from common import get_protein_name_dict, get_protein_coverage
import argparse
import bisect
from multiprocessing import Pool
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(
	description='Reads a peptide database, gives statistics on protein coverage by non-specific / protein-specific / proteoform-specific peptides.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input TSV file, peptide database")

parser.add_argument("-f", dest="fasta_file", required=True,
                    help="input FASTA file")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="# threads to use")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

print ('Reading', args.fasta_file)
name_dict, all_proteins = get_protein_name_dict(args.fasta_file)

print ('Reading', args.input_file)
pep_df = pd.read_csv(args.input_file, sep='\t', header=0)

result_columns = ['ProteinID', 'coverage']
result_data = []

# iterate through the peptide list - list is sorted by gene ID
current_gene = ''

total_aa = [0, 0, 0, 0]

all_peptides = []

# split up the DF by gene ID
all_genes = pep_df[['GeneID']].drop_duplicates()['GeneID'].tolist()

def get_gene_df(geneID):
    return pep_df[pep_df['GeneID'] == geneID]

#for geneID in all_genes:
with Pool(args.threads) as p:
    all_peptides = p.map(get_gene_df, all_genes)

def process_gene(idx):
    geneID = all_genes[idx]
    local_df = all_peptides[idx]
    
    protein_peptides = {}       # dictionary of mapptings between peptides and proteins -> to be aggregated into coverage stats
    local_aa = [0, 0, 0, 0]     # length of covered regions by type: 0: shoudl not be covered; 1: always canonical; 2: possibly protein-specific; 3: possibly proteoform-specific
    result = []

    # loop through peptides maping to this gene
    for index, row in local_df.iterrows():
        if geneID == '-' or row['matches_contaminant']:
            continue

        pep_length = len(row['Sequence'])

        # align this peptide to each matching protein (haplotype)
        for i, proteinID in enumerate(row['ProteinID'].split(',')):
            protein_name = name_dict[proteinID]
            ref_stable_id = protein_name.split(':')[0]

            pep_type = row['category']
            encoded_type = -1
            if pep_type == 'non-unique':
                encoded_type = 0
            elif pep_type == 'protein-specific':
                encoded_type = 1
            elif pep_type == 'proteoform-specific':
                encoded_type = 2

            pep_pos = int(row['Position'].split(',')[i])

            # add into a dictionary of peptide - ref. protein mappings
            if ref_stable_id in protein_peptides:
                protein_peptides[ref_stable_id]['peptides'].append([pep_pos, pep_length, encoded_type])
            else:
                protein_peptides[ref_stable_id] = { 'peptides': [ [pep_pos, pep_length, encoded_type] ] }

    # aggregate coverage by each protein
    for proteinID in protein_peptides:
        coverage = get_protein_coverage(sorted(protein_peptides[proteinID]['peptides'], key=lambda x: x[0]))

        for region in coverage:
            region_len = region[1] - region[0]
            region_type = region[2] + 1
            local_aa[region_type] += region_len

        coverage_str = ";".join(list(map(lambda x: "_".join(list(map(lambda y: str(y), x))), coverage)))

        result.append([proteinID, coverage_str])

    return [ result, local_aa ]

print ('Annotating proteome coverage.')

with Pool(args.threads) as p:
    perm = np.random.permutation(len(all_genes))
    gene_results = p.map(process_gene, perm)
    #gene_results = list(map(process_gene, all_genes))
    for i, result in enumerate(gene_results):
        if (len(result[0]) == 0):
            continue 

        result_data.extend(result[0])

        total_aa[0] += result[1][0]
        total_aa[1] += result[1][1]
        total_aa[2] += result[1][2]
        total_aa[3] += result[1][3]

print ('Done.')

result_df = pd.DataFrame(columns=result_columns, data=result_data)
result_df.to_csv(args.output_file, sep='\t', header=True, index=False)

#total_aa_sum = sum(total_aa)

# check the length of the proteome (= sum of lengths of all canonical proteins in fasta)
total_aa_sum = 0

for protein in all_proteins.values():
    if (protein['accession'].startswith('ENSP')):
        total_aa_sum += len(protein['sequence'])

print ("Proteome length:", total_aa_sum)
print ("Possible non-specific peptide: %d AAs - %.2f %%,\npossible protein-specific peptides: %d AAs - %.2f %%,\npossible proteoform-specific peptides: %d AAs - %.2f %%" % (total_aa[1], (total_aa[1] / total_aa_sum) * 100, total_aa[2], (total_aa[2] / total_aa_sum) * 100, total_aa[3], (total_aa[3] / total_aa_sum) * 100))
