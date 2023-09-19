import argparse
import re
import pandas as pd
import bisect
from multiprocessing import Pool
import numpy as np
from tqdm import tqdm
from common import read_fasta, digest

parser = argparse.ArgumentParser(
	description='Reads a Protein DB FASTA, creates a list of peptides for each gene with assicoated peptide types (reference/variant/haplotype)')

parser.add_argument("-i", dest="input_fasta", required=True,
                    help="input FASTA file")

parser.add_argument("-m", dest="missed_cl", type=int, required=True,
                    help="# missed cleavage sites allowed")

parser.add_argument("-min_len", dest="min_len", type=int, required=True,
                    help="min. peptide length")

parser.add_argument("-max_len", dest="max_len", type=int, required=True,
                    help="max. peptide length")

parser.add_argument("-sl", dest="subst_list", required=True,
                    help="AA substitution list by protein")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="# threads to use")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

print ("Creating a list of tryptic peptides from:", args.input_fasta, '\n# missed cleavage sites:', args.missed_cl)

all_proteins = read_fasta(args.input_fasta)
total_protein_count = len(all_proteins)
#proteins_processed = 0

subst_df = pd.read_csv(args.subst_list, sep='\t', converters={i: str for i in range(10)})
subst_df = subst_df.set_index('protein_stable_id')

peptides_columns = ['Sequence', 'GeneID', 'ProteinID', 'Position', 'possible_linked_SNPs', 'matches_contaminant', 'type', 'SNPs']
peptides_data = []
contaminant_peptides = []

possible_linked_snps = []

# helper for possible SNPs - how many SNPs fall into the peptide?
def get_ovelapping_snps(pep_position, pep_len, variant_starts, variant_ends):
        idx = bisect.bisect_left(variant_ends, pep_position)
        result_idxs = []

        for i in range(idx, len(variant_starts)):
                var_start = variant_starts[i]
                var_end = variant_ends[i]
                if (var_start > pep_position and var_start <= pep_position + pep_len):
                        result_idxs.append(i)
                elif (var_start <= pep_position and var_end > pep_position + 1):
                        result_idxs.append(i)
                #elif (var_start > pep_position + pep_len):
                #        break

        return result_idxs

# helper for possible SNPs - can the peptide contain a SNPs? (used for haplotypes without a multi-variant peptide)
def has_possible_snp(pep_position, pep_len, variant_starts, ref_allele_lengths):
        result = False

        for i in range(0, len(variant_starts)):
                var_start = int(variant_starts[i])
                var_end = var_start + int(ref_allele_lengths[i])
                if (var_start > pep_position and var_start <= pep_position + pep_len):
                        result = True
                        break
                elif (var_start <= pep_position and var_end > pep_position + 1):
                        result = True
                        break
                    
        return result

all_protein_ids = list(all_proteins.keys())

# digest all proteins into peptides and save
#for proteinID in all_proteins:
def process_protein(idx):
    proteinID = all_protein_ids[idx]
    sequence = all_proteins[proteinID]['sequence']
    peptides, peptide_positions = digest(sequence, args.missed_cl, args.min_len, args.max_len)
    result = []

    protein_changes = []
    change_locations = []
    var_alleles = []
    var_lengths = []
    var_ends = []
    refID = proteinID
    geneID = '-'
    if ('gene:' in all_proteins[proteinID]['description']):
        geneID = all_proteins[proteinID]['description'].split('gene:')[1].split('.')[0]

    # if there is any variation, remember where and what happens 
    if (proteinID.startswith('enshap')): 
        haplo_description = all_proteins[proteinID]['description'].split('haplotype:')[1]
        protein_changes = haplo_description.split(",")
        change_locations = list(map(lambda x: int(re.split('[a-zA-Z\*]+', x)[0]), protein_changes))        
        var_alleles = list(map(lambda x: re.split('\d+', x)[1].split('>')[1], protein_changes)) 

        # since no indels are included, alt and ref allele sequences have the same length
        var_lengths = list(map(lambda x: len(x), var_alleles))
        var_ends = [ change_locations[i] + var_len for i,var_len in enumerate(var_lengths) ]   

        refID = all_proteins[proteinID]['description'].split('ref:')[1].split('.')[0]

    # info about all possible mutations in this protein
    is_contaminant = not refID.startswith('ENSP')

    possible_mutations = None   # all mutations in all haplotypes of this protein

    if not is_contaminant:
            possible_mutations = subst_df.loc[refID]

    # assign metadata to peptides
    for pep_idx, pep in enumerate(peptides):
        pos = peptide_positions[pep_idx]
        pep_len = len(pep)
        found_changes = []

        # add contaminant peptides to list
        if is_contaminant and (pep not in contaminant_peptides):
            contaminant_peptides.append(pep)

        # could any SNPs happen in this region?
        possible_snp_count = []
        if (not is_contaminant) and (len(possible_mutations['subst_positions']) > 0):
            pos_lists = possible_mutations['subst_positions'].split(';')    # list of start positions of SNPS for all possible haplotypes 
            end_lists = possible_mutations['subst_ends'].split(';')         # list of end positions of SNPS for all possible haplotypes
            idxs_to_check = [ i for i,hap in enumerate(possible_mutations['possibly_haplotypic'].split(';')) if hap == '1' ]    # indexes of haplotypes that may protuce a multi-variant peptide

            # check only relevant haplotypes for multi-variant peptides
            for idx in idxs_to_check:
                pos_list = [ int (elem) for elem in pos_lists[idx].split(',') ]
                end_list = [ int (elem) for elem in end_lists[idx].split(',') ]

                zipped = list(zip(pos_list, end_list))
                zipped.sort(key=lambda x: x[1])
                pos_list, end_list = zip(*zipped)

                possible_snp_count.append(len(get_ovelapping_snps(pos, pep_len, pos_list, end_list)))

            # for the rest, just check if any of the SNPs could fall into the peptide
            possible_snp_count.append(int(has_possible_snp(pos, pep_len, possible_mutations['unique_snp_locations'].split(','), possible_mutations['unique_snp_lengths'].split(','))))
        else:
            possible_snp_count.append(0)

        # are any of the SNPs happening in this peptide?
        #   - sort changes in this haplotype by end position
        if (proteinID.startswith('enshap')): 
            zipped = list(zip(change_locations, var_ends, protein_changes))
            zipped.sort(key=lambda x: x[1])
            positions_sorted, ends_sorted, changes_sorted = zip(*zipped)

            #   - check which changes fall into this peptides (accounts for longer changes that start before the beginning of this peptide)
            found_idxs = get_ovelapping_snps(pos, pep_len, positions_sorted, ends_sorted)
            found_positions = [ positions_sorted[i] for i in found_idxs ]
            found_changes = [ changes_sorted[i] for i in found_idxs ]

            #   - sort back by start position
            if (len(found_changes) > 0):
                zipped = list(zip(found_positions, found_changes))
                zipped.sort(key=lambda x: x[0])
                found_positions, found_changes = zip(*zipped)

        pep_info = [pep, geneID, proteinID, pos, max(possible_snp_count), (is_contaminant or (pep in contaminant_peptides))]

        if (len(found_changes) == 1):
            pep_info.append('single_variant')
            pep_info.append(refID + ':' + found_changes[0])
        elif (len(found_changes) > 1):
            pep_info.append('multi_variant')
            pep_info.append(refID + ':' + ','.join(found_changes))
        else:
            pep_info.append('canonical')
            pep_info.append('-')

        result.append(pep_info)

    #proteins_processed += 1
    #print(str(proteins_processed) + ' / ' + str(total_protein_count), end='\r')
    return result

with Pool(args.threads) as p:
    #perm = np.random.permutation(total_protein_count)
    all_results = list(tqdm(p.imap_unordered(process_protein, range(0, total_protein_count)), total=total_protein_count))
    #all_results = list(map(process_protein, all_proteins.keys()))
    print ('Done, aggregating results.')
    for protein_peptides in all_results:
        for peptide_info in protein_peptides:
            peptides_data.append(peptide_info)
    del all_results

#all_results = list(map(process_protein, range(0, total_protein_count)))

peptides_df = pd.DataFrame(data=peptides_data, columns=peptides_columns)
peptides_df['Length'] = peptides_df['Sequence'].apply(lambda x: len(x))

del peptides_data

def group_rows(df):
    df['ProteinID'] = ','.join(df['ProteinID'].tolist())
    df['Position'] = df['Position'].apply(str)
    df['Position'] = ','.join(df['Position'].tolist())
    df['possible_linked_SNPs'] = ','.join(df['possible_linked_SNPs'].apply(str).tolist())
    df['matches_contaminant'] = any(df['matches_contaminant'].tolist())
    df['type'] = ','.join(df['type'].tolist())
    df['SNPs'] = '|'.join(df['SNPs'].tolist())

    return df

print ("Grouping by peptide and gene.")

def group_peptides(l):
    grouped_df = peptides_df[peptides_df['Length'] == l].groupby(['Sequence', 'GeneID']).apply(group_rows)
    grouped_df.drop_duplicates(subset=['Sequence', 'GeneID'], inplace=True)
    return grouped_df

#for l in range(args.min_len, args.max_len + 1):

with Pool(min(args.threads, (args.max_len - args.min_len + 1))):
    grouped_dfs = list(tqdm(p.imap_unordered(group_peptides, range(args.min_len, args.max_len + 1)), total=(args.max_len - args.min_len + 1)))
    print ("Writing output to", args.output_file)
    pd.concat(grouped_dfs).to_csv(args.output_file, sep='\t', header=True, index=False)
