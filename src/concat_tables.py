import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description='Concatinate all the provided tables into one TSV file, assuming all contain the same columns.')

parser.add_argument("-i", dest="input_filenames", required=True,
                    help="input files, comma-separated list")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output TSV file")

args = parser.parse_args()

inputs = args.input_filenames.split(",")
dataframes = []

for i,inputFilename in enumerate(inputs):
    df = pd.read_csv(inputFilename, sep='\t', header=0)
    dataframes.append(df)

result = pd.concat(dataframes)
result.to_csv(args.output_file, index=False, header=True, sep='\t')
