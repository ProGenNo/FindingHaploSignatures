from tqdm import tqdm
from common import get_protein_name_dict, get_protein_coverage
import argparse
import bisect
from multiprocessing import Pool
import pandas as pd

def parse_snp_list(SNPs):
    result = [] 
    for haplotypeGroup in SNPs.split(';'):
        for haplotypeName in haplotypeGroup.split('|'):
            if haplotypeName == '-':
                result.append('-')
                continue
            proteinID, SNPs = haplotypeName.split(':')
            for SNP in SNPs.split(','):
                result.append(proteinID + ':' + SNP)
    
    return result

def get_snp_stats(snp_list):
    counts = { 'single_variant': 0, 'both': 0, 'multi_variant': 0 }

    for snp in snp_list['multi_variant']:
        if snp in snp_list['single_variant']:
            counts['both'] += 1
        else:
            counts['multi_variant'] += 1

    counts['single_variant'] = len(snp_list['single_variant']) - counts['both']

    return counts

parser = argparse.ArgumentParser(
	description='Reads a peptide database, gives statistics on protein coverage by canonical/non-canonical peptides.')

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
#print ('Sorting peptides by gene ID.')
#pep_df.sort_values(by='GeneID', inplace=True)

total_peptide_count = len(pep_df)

result_columns = ['ProteinID', 'coverage', 'SNP_stats']
result_data = []

# iterate through the peptide list - list is sorted by gene ID
current_gene = ''
 # accessed by protein ID, 'peptides' list of lists, 3 values: start, length, type, 'var_covered': # of variants covered by variant / haplotypic peptides

total_aa = [0, 0, 0, 0]
total_var = [0, 0, 0]

all_peptides = {}

# split up the DF by gene ID
all_genes = pep_df[['GeneID']].drop_duplicates()['GeneID'].tolist()
for geneID in all_genes:
    all_peptides[geneID] = pep_df[pep_df['GeneID'] == geneID]

def process_gene(geneIdx):
    geneID = all_genes[geneIdx]
    local_df = all_peptides[geneID]
    
    protein_peptides = {}       # dictionary of mapptings between peptides and proteins -> to be aggregated into coverage stats
    local_aa = [0, 0, 0, 0]     # length of covered regions by type: 0: shoudl not be covered; 1: always canonical; 2: possibly single-variant; 3: possibly multi-variant
    local_var = [0, 0, 0]       # number of variants covered by peptide type: 0: only single-variant; 1: both; 2: only multi-variant
    result = []

    # loop through peptides maping to this gene
    for index, row in local_df.iterrows():
        if geneID == '-' or row['matches_contaminant']:
            continue

        pep_length = len(row['Sequence'])
        all_SNPs = parse_snp_list(row['SNPs'])

        # align this peptide to each matching protein (haplotype)
        for i, proteinID in enumerate(row['ProteinID'].split(',')):
            protein_name = name_dict[proteinID]
            ref_stable_id = protein_name.split(':')[0]

            pep_type = row['type'].split(',')[i]
            pep_type_global = row['type_aggregated']
            if pep_type_global != pep_type:
                continue

            encoded_type = -1
            if pep_type == 'canonical':
                encoded_type = 0
            elif pep_type == 'single_variant':
                encoded_type = 1
            elif pep_type == 'multi_variant':
                encoded_type = 2

            pep_pos = int(row['Position'].split(',')[i])

            SNPs = all_SNPs[i]
            if (SNPs != '-'):
                SNPs = SNPs.split(':')[1].split(',')
            else:
                SNPs = []

            # add into a dictionary of peptide - ref. protein mappings
            if ref_stable_id in protein_peptides:
                protein_peptides[ref_stable_id]['peptides'].append([pep_pos, pep_length, encoded_type])
            else:
                protein_peptides[ref_stable_id] = { 'peptides': [ [pep_pos, pep_length, encoded_type] ], 'var_covered': {'single_variant': [], 'multi_variant': []}  }

            # remember the variants found by peptide type
            if pep_type == 'multi_variant':
                for SNP in SNPs:
                    if SNP not in protein_peptides[ref_stable_id]['var_covered']['multi_variant']:
                        protein_peptides[ref_stable_id]['var_covered']['multi_variant'].append(SNP)
            elif pep_type == 'single_variant':
                for SNP in SNPs:
                    if SNP not in protein_peptides[ref_stable_id]['var_covered']['single_variant']:
                        protein_peptides[ref_stable_id]['var_covered']['single_variant'].append(SNP)

    # aggregate coverage by each protein
    for proteinID in protein_peptides:
        coverage = get_protein_coverage(sorted(protein_peptides[proteinID]['peptides'], key=lambda x: x[0]))
        SNPstats = get_snp_stats(protein_peptides[proteinID]['var_covered'])

        for region in coverage:
            region_len = region[1] - region[0]
            region_type = region[2] + 1
            local_aa[region_type] += region_len

        # add the remainder of the protein that's not covered, if any
        prot_len = len(all_proteins[proteinID]['sequence'])
        local_aa[0] += (prot_len - coverage[-1][1])
        coverage.append([coverage[-1][1], prot_len, -1])

        local_var[0] += SNPstats['single_variant']
        local_var[1] += SNPstats['both']
        local_var[2] += SNPstats['multi_variant']

        coverage_str = ";".join(list(map(lambda x: "_".join(list(map(lambda y: str(y), x))), coverage)))
        SNPstats_str = "single_variant:%d,both:%d,multi_variant:%d" % (SNPstats['single_variant'], SNPstats['both'], SNPstats['multi_variant'])

        result.append([proteinID, coverage_str, SNPstats_str])

    return [ result, local_aa, local_var ]

print ('Annotating proteome coverage.')

with Pool(args.threads) as p:
    gene_results = list(tqdm(p.imap_unordered(process_gene, range(len(all_genes))), total=len(all_genes)))
    p.close()
    p.join()

    for i, result in enumerate(gene_results):
        if (len(result[0]) == 0):
            continue 

        result_data.extend(result[0])

        total_aa[0] += result[1][0]
        total_aa[1] += result[1][1]
        total_aa[2] += result[1][2]
        total_aa[3] += result[1][3]
        
        total_var[0] += result[2][0]
        total_var[1] += result[2][1]
        total_var[2] += result[2][2]

print ('Done.')

result_df = pd.DataFrame(columns=result_columns, data=result_data)
result_df.to_csv(args.output_file, sep='\t', header=True, index=False)

#total_aa_sum = sum(total_aa)

# checl the length of the proteone (= sum of lengths of all canonical proteins in fasta)
total_aa_sum = 0

for protein in all_proteins.values():
    if (protein['accession'].startswith('ENSP')):
        total_aa_sum += len(protein['sequence'].replace('*', ''))

# total number of discoverable substitutions
total_var_sum = sum(total_var)

print ("Proteome length:", total_aa_sum)
print ("Canonical proteome: %d AAs - %.2f %%,\npossible single_variant peptides: %d AAs - %.2f %%,\npossible multi_variant peptides: %d AAs - %.2f %%,\nsequences not matching to peptides: %d AAs - %.2f %%" % (total_aa[1], (total_aa[1] / total_aa_sum) * 100, total_aa[2], (total_aa[2] / total_aa_sum) * 100, total_aa[3], (total_aa[3] / total_aa_sum) * 100, (total_aa_sum - sum(total_aa[1:])), ((total_aa_sum - sum(total_aa[1:])) / total_aa_sum) * 100))
print ()
print ("Total number of variants:", total_var_sum)
if (total_var_sum > 0):
    print ("Variants covered in var. pept.: %d - %.2f %%; \nvariants covered in both var. and hapl. pept.: %d - %.2f %%; \nvariants covered in multi_variant peptides: %d - %.2f %%" % (total_var[0], (total_var[0] / total_var_sum) * 100, total_var[1], (total_var[1] / total_var_sum) * 100, total_var[2], (total_var[2] / total_var_sum) * 100))
