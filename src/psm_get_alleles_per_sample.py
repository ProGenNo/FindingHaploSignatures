import argparse
import re
import pandas as pd

parser = argparse.ArgumentParser(
	description='Reads the PSM report with annotated variation, writes a new file with summarization of alleles per locus in the protein for each sample.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input TSV file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()

print ('Reading', args.input_file)
psm_df = pd.read_csv(args.input_file, sep='\t', header=0)
psm_df = psm_df[(psm_df['PeptideType'] != 'decoy') & (psm_df['PeptideType'] != 'contaminant') & (psm_df['PeptideType'] != 'has_stop')]

samples_data = {}

for index,row in psm_df.iterrows():
    if ((row['q-value'] > 1.01)): #or (('variant' in row['PeptideType']) and ('confident' not in row['PepQuery_class']))):
        continue

    sample = row['sample_ID']
    proteins = {}
    processed_snps = []

    if (row['CoveredSNPs'] != '-'):
        for substr in row['CoveredSNPs'].split(';'):
            protID = substr.split(':',1)[0]

            for SNP in substr.split(':',1)[1].split(','):
                if (SNP in processed_snps):
                    continue

                processed_snps.append(SNP)  # avoid duplicates
                loc = re.split('[^\d+]', SNP)[0]
                allele = SNP.split('>',1)[1]

                if (protID in proteins):
                    if (loc in proteins[protID]):
                        if allele in proteins[protID][loc]['alleles']:
                            allele_idx = proteins[protID][loc]['alleles'].index(allele)
                            proteins[protID][loc]['PSM_counts'][allele_idx] += 1
                        else:
                            proteins[protID][loc]['alleles'].append(allele)
                            proteins[protID][loc]['PSM_counts'].append(1)
                    else:
                        proteins[protID][loc] = {'alleles': [allele], 'PSM_counts': [1]}
                else:
                    proteins[protID] = {}
                    proteins[protID][loc] = {'alleles': [allele], 'PSM_counts': [1]}

    if (row['CoveredRefAlleles'] != '-'):
        for REF in re.split(r"[,;]", row['CoveredRefAlleles']):
            protID = REF.split(':',1)[0]
            loc = REF.split(':',2)[1]
            allele = REF.split(':',2)[2]

            if (protID in proteins):
                if (loc in proteins[protID]):
                    if allele in proteins[protID][loc]['alleles']:
                        allele_idx = proteins[protID][loc]['alleles'].index(allele)
                        proteins[protID][loc]['PSM_counts'][allele_idx] += 1
                    else:
                        proteins[protID][loc]['alleles'].append(allele)
                        proteins[protID][loc]['PSM_counts'].append(1)
                else:
                    proteins[protID][loc] = {'alleles': [allele], 'PSM_counts': [1]}
            else:
                proteins[protID] = {}
                proteins[protID][loc] = {'alleles': [allele], 'PSM_counts': [1]}

    if sample in samples_data:
        for protID in proteins:
            if (protID in samples_data[sample]):
                for loc in proteins[protID]:
                    if (loc in samples_data[sample][protID]):
                        for local_allele_idx,allele in enumerate(proteins[protID][loc]['alleles']):
                            if allele in samples_data[sample][protID][loc]['alleles']:
                                allele_idx = samples_data[sample][protID][loc]['alleles'].index(allele)
                                samples_data[sample][protID][loc]['PSM_counts'][allele_idx] += proteins[protID][loc]['PSM_counts'][local_allele_idx]
                            else:
                                samples_data[sample][protID][loc]['alleles'].append(allele)
                                samples_data[sample][protID][loc]['PSM_counts'].append(proteins[protID][loc]['PSM_counts'][local_allele_idx])
                    else:
                        samples_data[sample][protID][loc] = proteins[protID][loc]
            else:
                samples_data[sample][protID] = proteins[protID]
    else:
        samples_data[sample] = proteins

result_data = []
result_columns = [
    'sample',
    'proteinID',
    'location',
    'alleles',
    'psms_per_allele'
]

for sampleID in samples_data:
    for protID in samples_data[sampleID]:
        for loc in samples_data[sampleID][protID]:
            result_data.append([sampleID, protID, loc, ';'.join(samples_data[sampleID][protID][loc]['alleles']), ';'.join([ str(i) for i in samples_data[sampleID][protID][loc]['PSM_counts']])])

result_df = pd.DataFrame(data=result_data, columns=result_columns)
result_df.to_csv(args.output_file, sep='\t', index=False, header=True)


