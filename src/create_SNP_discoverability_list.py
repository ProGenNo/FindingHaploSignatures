import argparse
import bisect
import pandas as pd


parser = argparse.ArgumentParser(
        description='Reads the peptide list TSV file with haplotype frequency column, creates a list of SNPs and their discoverability by peptide class and highest haplotype frequency.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input tab-separated file")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

args = parser.parse_args()


pep_df = pd.read_csv(args.input_file, sep='\t', header=0)

# three lists sorted in the same way
hashes = []
snps = []
pep_types = []
freq_differences = []

for i in range(len(pep_df)):
	row = pep_df.iloc[i]
	#pep_type = row['type']
	freq_diff = -1
	if row['highest_frequency'] == row['highest_frequency']:				# check for NaN values!
		hap_freq = float(row['highest_frequency'].split(':', 1)[1].split('(', 1)[0])
		ref_freq = float(row['highest_frequency'].split('(', 1)[1].split(')', 1)[0])
		freq_diff = hap_freq - ref_freq
	else:
		continue

	for i,SNPlist in enumerate(row['SNPs'].split(';')):
		typelist = row['type'].split(';')[i].split(',')

		for j,SNPgroup in enumerate(SNPlist.split('|')):
			pep_type = typelist[j]
			refID = SNPgroup.split(':', 1)[0]

			for SNP in SNPgroup.split(':', 1)[1].split(','):
				SNPID = refID + ':' + SNP
				h = hash(SNPID)

				nearest_idx = bisect.bisect_left(hashes, h)
				if (nearest_idx >= len(hashes) or hashes[nearest_idx] != h):
					hashes.insert(nearest_idx, h)
					snps.insert(nearest_idx, SNPID)
					pep_types.insert(nearest_idx, [pep_type])
					freq_differences.insert(nearest_idx, freq_diff)
				else:
					if pep_type not in pep_types[nearest_idx]:
						pep_types[nearest_idx].append(pep_type)

					if freq_diff > freq_differences[nearest_idx]:
						freq_differences[nearest_idx] = freq_diff

results_df = pd.DataFrame()
results_df['SNP'] = snps
results_df['pep_types'] = pep_types
results_df['freq_diffs'] = freq_differences

results_df.to_csv(args.output_file, sep='\t', header=True, index=False)
