import pandas as pd
import spectrum_utils.spectrum as sus
from pyteomics import mzml
import argparse
import re
import math

def parse_annotation(elem):
    if elem is None:
        return "None"

    return elem.annotation

def parse_annotation_pred(row):
    return row['ion'].lower() + str(row['ionnumber']) + '+'

parser = argparse.ArgumentParser(
	description='Export observed and predicted spectra.')

parser.add_argument("-i", dest="input_file", required=True,
                    help="input tab-separated file")

parser.add_argument("-m", dest="mzml_dir", required=True,
                    help="directory with mzml files")

parser.add_argument("-p", dest="predictions_dir", required=True,
                    help="directory with prediction files")

parser.add_argument("-e", dest="engine", required=False,
                    help="search engine name (default: \"xtandem\")", default='xtandem')

#parser.add_argument("-type", dest="pep_type", required=True,
#                    help="peptide type")

parser.add_argument("-o", dest="output_file", required=True,
                    help="output filename")

args = parser.parse_args()

psm_df = pd.read_csv(args.input_file, sep='\t', header=0)
#psm_df = psm_df[psm_df['PeptideType'] == args.pep_type]

mzml_readers = {}

psm_df.sort_values(by=['posterior_error_prob'], inplace=True)

fragment_tol_mass = 10
fragment_tol_mode = 'ppm'

columns = ['PSMId', 'peptide', 'spectrum_title', 'posterior_error_prob','observed_mz', 'observed_intensity', 'observed_annotation', 'predicted_mz', 'predicted_intensity', 'predicted_annotation']
data = []

i = 0

for index, row in psm_df.iterrows():

    # get the identification charge from the 1 hot encoding
    ident_charge = -1 
    if row['charge_2']:
        ident_charge = 2
    elif row['charge_3']:
        ident_charge = 3
    elif row['charge_4']:
        ident_charge = 4

    # check for missing info - skip if missing
    if ( math.isnan(row['measured_mz']) or (ident_charge == -1)):
        print ('Skipping row', i, '(missing data)')
        i = i + 1
        continue

    sample = row['fragment_ID']
    if sample not in mzml_readers:
        mzml_readers[sample] = mzml.read(args.mzml_dir + sample + '.mzML')
    
    spectrum_title = row['SpectrumTitle']
    spectrum_pred_id = int(row['PSMId'].split('_')[1])
    peptide = row['Sequence']

    # clean up peptide sequence
    #print('Peptide:', peptide)
    pattern = r'[\W\d]'
    peptide = re.sub(pattern, '', peptide)
    #print('Peptide:', peptide)


    precursor_mass = row['measured_mz']
    
    print (i, '/', len(psm_df), 'Processing:', sample, 'peptide', peptide, 'Spectrum prediction id', spectrum_pred_id)
    print ('Identification Charge', ident_charge, 'Theoretical Mass', precursor_mass)

    # ids_df = pd.read_csv(args.id_dir + sample + '/' + sample + "_" + args.engine + "_PSM_identifiers_export", sep='\t', header=0)
    ms2pipId = spectrum_pred_id

    pred_df = pd.read_csv(args.predictions_dir + sample + "_" + args.engine + "/peprec_HCDch2_predictions.csv", header=0)

    pred_df = pred_df[pred_df['spec_id'] == ms2pipId]
    pred_df = pred_df[pred_df['prediction'] > 0.001]
    #print(pred_df)
    pred_df['annotation'] = pred_df.apply(parse_annotation_pred, axis=1)
    #print(pred_df['annotation'])

    spectrum = mzml_readers[sample].get_by_id(spectrum_title)

    spectrum_top = sus.MsmsSpectrum("my_spectrum", precursor_mass, ident_charge, spectrum['m/z array'], spectrum['intensity array'], peptide=peptide)
    #spectrum_bottom = sus.MsmsSpectrum("my_prediction", precursor_mass, ident_charge, pred_df['mz'].tolist(), pred_df['prediction'].tolist(), peptide=peptide)

    min_mz = min([ min(spectrum['m/z array']), min(pred_df['mz'].tolist()), 100 ])
    max_mz = max([ max(spectrum['m/z array']), max(pred_df['mz'].tolist()), 1400 ])

    print('Max. m/z:', max_mz, 'Min. m/z:', min_mz)

    spectrum_top = (spectrum_top.set_mz_range(min_mz=min_mz, max_mz=max_mz)
                    .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
                    .filter_intensity(min_intensity=0.05, max_num_peaks=50)
                    .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
                                    ion_types='by'))

    #spectrum_bottom = (spectrum_bottom.set_mz_range(min_mz=min_mz, max_mz=max_mz)
    #                .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
    #                .filter_intensity(min_intensity=0.05, max_num_peaks=50)
    #                .annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
    #                                ion_types='by'))

    data.append([row['PSMId'], peptide, spectrum_title, row['posterior_error_prob'], spectrum_top.mz.tolist(), spectrum_top.intensity.tolist(), list(map(parse_annotation, spectrum_top.annotation.tolist())), pred_df['mz'].tolist(),  pred_df['prediction'].tolist(), pred_df['annotation'].tolist() ])
    i = i + 1

values_df = pd.DataFrame(data=data, columns=columns)

print ('Writing output file:', args.output_file)
values_df.to_csv(args.output_file, sep='\t', header=True, index=False)
