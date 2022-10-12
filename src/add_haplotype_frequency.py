import pandas as pd
import numpy as np
import argparse
from multiprocessing import Pool

parser = argparse.ArgumentParser(
	description='Adds the FoO of the most frequent haplotype matching this peptide, and the FoO of the reference haplotype.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input TSV file")

parser.add_argument("-hap", dest="haplotypes_file", required=True,
                    help="protein haplotypes CSV file")

parser.add_argument("-t", dest="threads", type=int, required=True,
                    help="# threads to use")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output file")

parser.add_argument("-var", dest="var_list_file", required=False,
                    help="optional: file to list the discoverable AA variants")

args = parser.parse_args()

df = pd.read_csv(args.input_file, sep='\t', header=0)
df = df[~(df['type'].str.contains('reference'))]

haplo_df = pd.read_csv(args.haplotypes_file, header=0)

# get the frequency of the most frequent haplotype that contains this peptide
# and compare with the frequency of the reference
def get_freq(idx):
    row = df.iloc[idx]
    freqs = []

    for haplotypeGroup in row['SNPs'].split(';'):
        ref_freq = -1
        matching_hapl_freq = -1

        for haplotypeName in haplotypeGroup.split('|'):
            proteinID, SNPs = haplotypeName.split(':')
            
            relevantHaplotypes = haplo_df[haplo_df["protein_stable_id"].str.match(proteinID)]
            refHaplotype = relevantHaplotypes[relevantHaplotypes["name"].str.match(proteinID + ":REF")]
            try:
                matchingHaplotypes = relevantHaplotypes[relevantHaplotypes["name"].str.contains(proteinID) & relevantHaplotypes["name"].str.contains(SNPs)]
            except:
                print('ProteinID:', proteinID, 'SNPs:', SNPs)
                return ""

            if (len(refHaplotype) > 0):
                ref_freq = refHaplotype[["freq"]].values.tolist()[0][0]
            if (len(matchingHaplotypes) > 0):
                matching_hapl_freq = max(list(map(lambda x: x[0], matchingHaplotypes[["freq"]].values.tolist())))
	
        if (matching_hapl_freq <= ref_freq):
            freqs.append("less:" + str(matching_hapl_freq) + "(" + str(ref_freq) + ")")
        else:
            freqs.append("more:" + str(matching_hapl_freq) + "(" + str(ref_freq) + ")")

    return ";".join(freqs)

#df['Highest Frequency'] = df.apply(get_freq, axis=1)
with Pool(args.threads) as p:
    perm = list(np.random.permutation(len(df)))
    results = p.map(get_freq, perm)

    # sort back the results
    zipped = list(zip(perm, results))
    zipped.sort(key=lambda x: x[0])
    perm, results = zip(*zipped)

    df['highest_frequency'] = results
    df.to_csv(args.output_file, sep='\t', header=True, index=False)

# Create a list of AA variants that can be found in haplotypic peptides, and assign the frequency of the haplotype compared to reference
if (args.var_list_file):
    vars_data = []

    for index, row in df.iterrows():
        for haplotypeGroup in row['SNPs'].split(';'):
            for haplotypeName in haplotypeGroup.split('|'):
                proteinID, SNPs = haplotypeName.split(':')
                for SNP in SNPs.split(','):
                    vars_data.append([proteinID + ':' + SNP, row['highest_frequency']])

    df_vars = pd.DataFrame(data=vars_data, columns=["variant", "highest_haplotype_frequency"])
    df_vars.to_csv(args.var_list_file, sep='\t', header=True, index=False)
