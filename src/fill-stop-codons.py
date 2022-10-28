import argparse
from common import read_fasta

parser = argparse.ArgumentParser(description='Reads a FASTA file, fills in missing stop codons from the info in the header')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input FASTA file", metavar="FILE")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output FASTA file", metavar="FILE")

args = parser.parse_args()

all_proteins = read_fasta(args.input_file)

outfile = open(args.output_file, 'w')

for protein in all_proteins.values():
    sequence = protein['sequence']
    if 'stop:' in protein['description']:
        positions = [ int(p) for p in protein['description'].split('stop:', 1)[1].split(maxsplit=1)[0].split('.') ]
        for pos in positions:
            sequence = sequence[:pos] + '*' + sequence[pos:]

    outfile.write('>' + protein['tag'] + '|' + protein['accession'] + '|' + protein['description'] + '\n')
    outfile.write(sequence + '\n')

outfile.close()
