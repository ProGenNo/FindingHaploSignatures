import argparse
import re
import pandas as pd
import bisect
from tqdm import tqdm
from multiprocessing import Pool
from common import get_protein_name_dict

parser = argparse.ArgumentParser(
	description='Reads the PSM report file, creates a new file with information about covered SNPs for each PSM and protein')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input tab-separated file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

parser.add_argument("-hap", dest="haplo_db", required=True,
                    help="haplotypes csv file")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="maximum number of threads")

parser.add_argument('-s', dest="sample", required=True,
                    help="Sample ID, or 'all' if processing a merged PSM report")  

parser.add_argument("-f", dest="fasta_file", required=True,
                    help="fasta file")                    

args = parser.parse_args()

print ("Reading", args.fasta_file)
protein_name_dict, all_proteins = get_protein_name_dict(args.fasta_file)

print ("Reading", args.input_file)
psm_df = pd.read_csv(args.input_file, sep='\t', header=0)
psm_count = len(psm_df)

print ("Reading", args.haplo_db)
haplo_df = pd.read_csv(args.haplo_db, header=0)

# Aggregate all the reference allele locations for each protein
ref_alleles = {}
for index,row in haplo_df.iterrows():
    if (('REF' in row['name']) or ('ins' in row['name']) or ('del' in row['name']) or ('{' in row['name'])):
        continue

    protID = row['protein_stable_id']
    SNPs = row['name'].split(':',1)[1].split(',')
    try:
        locations = [ int(re.split('[^\d]+', SNP)[0]) for SNP in SNPs ]
        ref = [ re.split('\d+', SNP)[1].split('>',1)[0] for SNP in SNPs ]
    except:
        print(row['name'])
    #ref = [ re.split('\d+', SNP)[1].split('>',1)[0] for SNP in SNPs ]
    
    if protID in ref_alleles:
        protData = ref_alleles[protID]

        for i in range(len(SNPs)):
            idx = bisect.bisect_left(protData['loc'], locations[i])
            if ((idx == len(protData['loc'])) or (locations[i] != protData['loc'][idx])):
                protData['loc'].insert(idx, locations[i])
                protData['ref'].insert(idx, ref[i])

        ref_alleles[protID] = protData
    
    else:
        ref_alleles[protID] = { 'loc': locations, 'ref': ref }

annotation_data = []
annotation_columns = [
        'PSMId', 
        'Proteins',     
        'Position',     # Replace the positions of peptides within proteins taking into account possible preceding stop codons, which have been removed in the search database, but need to be included in the post-processing step to correctly align sequences
        'PeptideType', 
        'CoveredSNPs', 
        'CoveredRefAlleles', 
        'Haplotypes', 
        'HaplotypeFreqs', 
        'OtherMatches',
        'GeneIDs'
        ]

print ("Annotating PSMs:")

# check if a given region of the protein contains a stop codon
# return boolean value, and position of the region within the whole protein when preceding stop codons are included
def check_stop(fastaID, pos_from, pos_to):
    desc = all_proteins[fastaID.split('.')[0]]['description']

    if 'stop:' not in desc:
        return False, pos_from

    stop_positions = sorted([ int(sp) for sp in desc.split('stop:',1)[1].split('.') ])

    for sp in stop_positions:
        if sp < pos_from:   # stop is before the region
            pos_from += 1
            pos_to += 1
        elif sp < pos_to:   # stop is in the region
            return True, pos_from
        else:
            break

    return False, pos_from

# Checks if any reference alleles are identified for the given protein
def check_ref_alleles(seq, protID, pos_from, pos_to):
    if (protID not in ref_alleles):
        return []
    result = []
    protData = ref_alleles[protID]
    for i,loc in enumerate(protData['loc']):
        if (loc >= pos_to):
            break
        if (loc >= pos_from):
            if (seq[loc - pos_from - 1] == protData['ref'][i]):
                result.append(protID + ':' + str(loc) + ':' + protData['ref'][i])

    return result

def process_row(index):
    row = psm_df.iloc[index]
    #print(row['Position'])

    fasta_peptide_starts = re.split(r"[,;]", row['Position'])
    peptide_starts = []
    protein_accessions = []

    for i,acc in enumerate(re.split(r"[,;]", row['Proteins'])):
        positions = [ int(pos) for pos in fasta_peptide_starts[i].split(',')]
        for pos in positions:
            protein_accessions.append(acc)
            peptide_starts.append(pos)

    #protein_accessions = row['Proteins'].split(';')      # all the accessions of candidate proteins
    peptide_length = len(row['Sequence'])                # length of the peptide sequence
    #peptide_starts = list(map(lambda x: int(x), row['Position'].replace(',', ';').split(';')))    # positions of the peptide within respective candidate protein

    is_decoy = all([ acc.endswith('REVERSED') for acc in protein_accessions ])
    peptide_starts = [ peptide_starts[i] for i,acc in enumerate(protein_accessions) if not acc.endswith('REVERSED') ]
    protein_accessions = [ acc for acc in protein_accessions if not acc.endswith('REVERSED') ]

    # check if peptide matches to a contaminant - discard if so
    is_contaminant = not all(list(map(lambda x: x.startswith('enshap') or x.startswith('ensvar') or x.startswith('ENSP'), protein_accessions)))

    # filter out matches that include a stop codon
    # adjust the peptide starts so that any preceding stop codons are taken into account
    have_stop = []
    for i,fastaID in enumerate(protein_accessions):
        has_stop, adjusted_start = check_stop(fastaID, peptide_starts[i], peptide_starts[i] + peptide_length)
        have_stop.append(has_stop)
        peptide_starts[i] = adjusted_start

    peptide_starts = [ peptide_starts[i] for i,has_stop in enumerate(have_stop) if not has_stop ]
    protein_accessions = [ protein_accessions[i] for i,has_stop in enumerate(have_stop) if not has_stop ]

    all_have_stop = all(have_stop)

    peptide_types = []
    protein_haplotypes = []
    protein_variants = []
    other_proteins = []
    gene_IDs = [ all_proteins[acc]['description'].split('gene:')[1].split('.')[0] for acc in protein_accessions if ('gene:' in all_proteins[acc]['description']) ]
    SNPs = []
    all_ref_alleles = []
    haplotype_frequencies = []
    pep_type = '-' # here will be specified whether the peptide is haplotype- or variant-specific

    # check if the peptide matches to a canonical protein - if so, do not annotate haplotypes
    is_canonical = any(list(map(lambda x: x.startswith('ENSP'), protein_accessions)))

    if (not is_canonical) and (not is_contaminant) and (not is_decoy) and (not all_have_stop):
        all_matches = []

        for prot_idx, acc in enumerate(protein_accessions):
            if ('enshap' not in acc or 'REVERSED' in acc or acc not in all_proteins):
                peptide_types.append("-")
                other_proteins.append(acc)
                all_matches.append(acc)
                continue

            haplo_description = all_proteins[acc]['description'].split('haplotype:')[1]
            protein_id = all_proteins[acc]['description'].split('ref:')[1].split('.')[0]
            protein_changes = haplo_description.split(",")
            haplotype_name = protein_name_dict[acc]
            all_matches.append(haplotype_name)
            all_ref_alleles.extend(check_ref_alleles(row['Sequence'], protein_id, peptide_starts[prot_idx], peptide_starts[prot_idx] + peptide_length))

            # get the reference and haplotype frequencies
            relevantHaplotypes = haplo_df[haplo_df["protein_stable_id"] == protein_id]
            refHaplotype = relevantHaplotypes[relevantHaplotypes["name"] == (protein_id + ":REF")]
            foundHaplotype = relevantHaplotypes[relevantHaplotypes["name"] == haplotype_name]

            haplo_freq = -1
            ref_freq = -1
            if (len(foundHaplotype) > 1):
                print('Haplotype name not unique, cannot get frequency:', haplotype_name)
            elif (len(foundHaplotype) == 0):
                print('Haplotype', haplotype_name, 'not found in database.')
            else:
                haplo_freq = foundHaplotype[["freq"]].values.tolist()[0][0]

            if (len(refHaplotype) > 0):
                ref_freq = refHaplotype[["freq"]].values.tolist()[0][0]

            freq_annot = ''
            if (haplo_freq <= ref_freq):
                freq_annot = "less:" + str(haplo_freq) + "(" + str(ref_freq) + ")"
            else:
                freq_annot = "more:" + str(haplo_freq) + "(" + str(ref_freq) + ")"

            # parse the haplotype name into AAs in reference and alt allele and locations of variation
            # ref_alleles = list(map(lambda x: re.split('\d+', x)[1].split('>')[0], protein_changes))
            # var_alleles = list(map(lambda x: re.split('\d+', x)[1].split('>')[1], protein_changes))
            change_locations = list(map(lambda x: int(re.split('[a-zA-Z\*]+', x)[0]), protein_changes))

            # check how many of these changes affect this peptide
            #found_changes_count = 0
            found_changes = []

            for change_idx, loc in enumerate(change_locations):
                try:
                    if (loc > peptide_starts[prot_idx] and loc <= peptide_starts[prot_idx] + peptide_length):
                        #print (acc + ":", "expected:", protein_changes[change_idx], "found:", row['Sequence'][:loc - peptide_starts[prot_idx - 1]], row['Sequence'][loc - peptide_starts[prot_idx] - 1], row['Sequence'][loc - peptide_starts[prot_idx]:])
                        #found_changes_count += 1
                        found_changes.append(protein_changes[change_idx])
                except:
                    raise ValueError(row['PSMId'], acc, "# peptide start entries:", len(peptide_starts), peptide_starts, '# protein accessions:', len(protein_accessions))

            if (len(found_changes) == 1):
                SNPs.append(protein_id + ':' + found_changes[0])
                protein_haplotypes.append(haplotype_name)
                haplotype_frequencies.append(freq_annot)

                if (pep_type != 'canonical'):
                    pep_type = 'single-variant'

            elif (len(found_changes) > 1):
                SNPs.append(protein_id + ':' + ','.join(found_changes))
                protein_haplotypes.append(haplotype_name)
                haplotype_frequencies.append(freq_annot)

                if (pep_type != 'canonical' and pep_type != 'single-variant'):
                    pep_type = 'multi-variant'
            else:
                other_proteins.append(haplotype_name)
                pep_type = 'canonical'

    else:
        # if the peptide matches to a canonical protein, discard all the haplotype matches
        other_proteins = [acc for acc in protein_accessions if acc.startswith("ENSP")]
        all_matches = protein_accessions        

        if (is_canonical):
            for i,protID in enumerate(protein_accessions):
                if protID.startswith("ENSP"):
                    all_ref_alleles.extend(check_ref_alleles(row['Sequence'], protID.split('.',1)[0], peptide_starts[i], peptide_starts[i] + peptide_length))

        if is_decoy:
            pep_type = 'decoy'
        elif all_have_stop:
            pep_type = 'has_stop'
        elif is_contaminant:
            pep_type = 'contaminant'
        else:
            pep_type = 'canonical'

    #psm_df.at[index, 'Peptide type'] = ','.join(peptide_types)
    if (len(SNPs) == 0):
        SNPs = ['-']
    if (len(all_ref_alleles) == 0):
        all_ref_alleles = ['-']
    if (len(protein_haplotypes) == 0):
        protein_haplotypes = ['-']
    if (len(haplotype_frequencies) == 0):
        haplotype_frequencies = ['-']
    if (len(other_proteins) == 0):
        other_proteins = ['-']
    if (len(gene_IDs) == 0):
        gene_IDs = ['-']

    #print(str(index) + ' / ' + str(psm_count), end='\r')
    return [row['PSMId'], ';'.join(protein_accessions), ';'.join([ str(pos) for pos in peptide_starts ]), pep_type, ';'.join(SNPs), ';'.join(all_ref_alleles),';'.join(protein_haplotypes), ';'.join(haplotype_frequencies), ';'.join(other_proteins), ';'.join(gene_IDs)]
    

# read PSMs line by line, check whether any peptide is haplotype- or variant-specific
# for index, row in psm_df.iterrows():
with Pool(args.threads) as pool:
    annotation_data = list(tqdm(pool.imap_unordered(process_row, range(0, len(psm_df))), total=len(psm_df)))
    pool.close()
    pool.join()

#annotation_data = list(map(process_row, range(0, len(psm_df)))) # sequential run for debugging

# there will be some empty rows - filter them out
annotation_data = [ row for row in annotation_data if len(row) > 0]

print ('Done')
#psm_df.to_csv("tonsil_trypsin_PSM_report_checked.txt", sep='\t', header=True, index=False)

annotation_df = pd.DataFrame(data=annotation_data, columns=annotation_columns)
result_df = annotation_df.join(psm_df.drop(['Proteins', 'Position'], axis=1).set_index('PSMId'), on='PSMId', rsuffix='_duplicate')
#summary_df.sort_values(by=['Peptide'], inplace=True)

print ('Writing results to', args.output_file)
result_df.to_csv(args.output_file, sep='\t', header=True, index=False)
