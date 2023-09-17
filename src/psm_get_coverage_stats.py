from common import read_fasta
import argparse
import pandas as pd
import bisect

# hepler class for bisect to compare objects
class KeyWrapper:
    def __init__(self, iterable, key):
        self.it = iterable
        self.key = key

    def __getitem__(self, i):
        return self.key(self.it[i])

    def __len__(self):
        return len(self.it)

    def insert(self, index, item):
        self.it.insert(index, item)

# map a region from the PSM coverage report to the predicted coverage by tryptic peptides
# returns length of coverage in region for classes: 0: should not be covered; 1: always canonical; 2: possibly single-variant; 3: possibly multi-variant
def get_region_class(start, end, all_regions):
    region_idx = max(bisect.bisect_left(KeyWrapper(all_regions, key=lambda x: x[0]), start) - 1, 0)
    region = all_regions[region_idx]
    result = [0, 0, 0, 0]

    if (start >= region[0]) and (end <= region[1]):
        result[region[2]+1] += end - start

    # peptide is overlapping a region boundary -> split
    elif (start >= region[0]):
        start_tmp = start
        result[region[2]+1] += region[1] - start_tmp
        start_tmp = region[1]

        while (region_idx < len(all_regions)):
            region = all_regions[region_idx]

            if (end < region[1]):
                result[region[2]+1] += end - start_tmp
                break

            result[region[2]+1] += region[1] - start_tmp
            start_tmp = region[1]
            region_idx += 1

    # print ('Peptide exceeding region boundaries.')
    return result

# avoid overlapping regions -> aggregate
def aggregate_regions(region_list):
    eventQ = []
    result = []

    for i,region in enumerate(region_list):
        eventQ.append({ 'type': 'start', 'pos': region[0]})
        eventQ.append({ 'type': 'end', 'pos': region[1]})

    eventQ.sort(key=lambda x: x['pos'])

    current_start = eventQ[0]['pos']
    active_regions = 1

    for i in range(1, len(eventQ)):
        event = eventQ[i]

        if event['type'] == 'start':
            if (active_regions == 0):
                current_start = event['pos']
            active_regions += 1

        else:
            active_regions -= 1
            if (active_regions == 0):
                result.append([current_start, event['pos']])
            elif (active_regions < 0):
                print ('Negative number of active regions')

    return result


parser = argparse.ArgumentParser(
	description='Reads a precomputed peptide stats file and the proteome coverage in PSMs, gives summary statistics on protein coverage by canonical/non-canonical peptides compared to coverage by discoverable peptides.')

parser.add_argument("-pep", dest="pep_stats", required=True,
                    help="input TSV file, peptide stats by protein")

parser.add_argument("-psm", dest="psm_coverage", required=True,
                    help="input TSV file, proteome coverage by PSMs")

parser.add_argument("-f", dest="fasta_file", required=True,
                    help="input FASTA file")

args = parser.parse_args()

all_proteins = read_fasta(args.fasta_file)

pep_df = pd.read_csv(args.pep_stats, sep='\t', header=0)
coverage_df = pd.read_csv(args.psm_coverage, sep='/t', header=0)

total_aa = [0, 0, 0, 0] # lengths of proteome covered by type

protein_all_regions = None        # list of tuples describing regions of the current protein [start, end, category]
current_protein_id = ""
covered_regions = []

for index, row in coverage_df.iterrows():
    proteinID = row['ProteinID']

    # after all regions of a protein (incl. its haplotypes) have been read -> aggregate and analyze
    if (current_protein_id != proteinID):
        if (current_protein_id != ""):
            aggregated_regions = aggregate_regions(covered_regions)

            for region in aggregated_regions:
                coverage = get_region_class(region[0], region[1], protein_all_regions)
                #if (coverage[0] > 0):
                #    print ('Protein', current_protein_id, 'region:', region, 'should not have been covered')

                #total_aa[0] += coverage[0]
                total_aa[1] += coverage[1]
                total_aa[2] += coverage[2]
                total_aa[3] += coverage[3]

        # cache predicted coverage for a new protein if needed
        protein_df = pep_df[pep_df['ProteinID'] == proteinID]
        if len(protein_df) < 1:
            print ('Protein', proteinID, 'not in database!')
            protein_all_regions = []
            current_protein_id = ""
        else:
            protein_all_regions = [ [ int(x) for x in region.split('_') ] for region in protein_df.iloc[0]['coverage'].split(';') ]
            current_protein_id = proteinID
        covered_regions = []

    covered_regions.extend([ [ int(x.split('(')[0]) for x in region.split('-') ] for region in row['CoveredRegions'].split(',') ])

# process the last protein
if (current_protein_id != ""):
    aggregated_regions = aggregate_regions(covered_regions)

    for region in aggregated_regions:
        coverage = get_region_class(region[0], region[1], protein_all_regions)
        #if (coverage[0] > 0):
        #    print ('Protein', current_protein_id, 'region:', region, 'should not have been covered')

        #total_aa[0] += coverage[0]
        total_aa[1] += coverage[1]
        total_aa[2] += coverage[2]
        total_aa[3] += coverage[3]

total_aa_sum = 0

for protein in all_proteins.values():
    if (protein['accession'].startswith('ENSP')):
        total_aa_sum += len(protein['sequence'])

print ("Proteome length:", total_aa_sum)
print ("Canonical proteome: %d AAs - %.2f %%,\npossible single-variant peptides: %d AAs - %.2f %%,\npossible multi-variant peptides: %d AAs - %.2f %%" % (total_aa[1], (total_aa[1] / total_aa_sum) * 100, total_aa[2], (total_aa[2] / total_aa_sum) * 100, total_aa[3], (total_aa[3] / total_aa_sum) * 100))
print ()
