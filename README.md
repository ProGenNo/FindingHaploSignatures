# Identifying Protein Haplotypes by Mass Spectrometry
Code related to the "Identifying Protein Haplotypes by Mass Spectrometry" publication: https://doi.org/10.1101/2022.11.21.517096

## Requirements and Usage
Required software is Snakemake and Conda, remaining libraries are included in the provided Conda environment, created automatically by Snakemake.

Steps for reproducing results:
- Download supplementary data from https://doi.org/10.6084/m9.figshare.21408117.v1
- Clone this repository
- Provide path to the downloaded files following instructions `config.yaml`
- For generating mirrored plots, a rerun of MS2PIP is necessary. This step can be skipped to re-generate the remaining results.
- Run the pipeline using `snakemake -c<# cores> -p --use-conda`

## Additional files provided
Two additional files were derived from data made public by Spooner et al. \[1\]

- `protein_haplotypes_1.csv`: List of protein haplotypes created using alleles with minor allele frequency (MAF) > 0.01
- `substitutions_list.tsv`: List of substitutions found in the provided haplotypes, along with other pre-computed values 

## Supplementary data format
Description of supplementary data accessible at https://doi.org/10.6084/m9.figshare.21408117.v1:

SD1 and SD2: protein databases in the FASTA format:
```
>tag|accession|description
PROTEINSEQUENCE
```
SD3: CSV, list of all resulting PSMs with confidence metrics. The included columns are:
- PSMId: artificial identifier of the peptide-spectrum match
- SpectrumTitle: identifier of the matching spectrum, as given in the data set
- Proteins: accession identifiers of the matching protein sequences
- Position: positions of the peptide within the respective protein
- Sequence: sequence of the peptide without post-translational modifications (PTMs)
- Peptide: sequence of the peptide with possible PTMs
- sample_ID: identifier of the sample as given in the data set
- fragment_ID: identifier of the sample fragment as given in the data set
- measured_rt: measured retention time (RT)
- rt_Abs_error: absolute distance between measured and predicted RT
- spectra_cos_similarity: cosine similarity between matching peaks in the observed and predicted spectrum
- spectra_angular_similarity: angular similarity between matching peaks in the observed and predicted spectrum
- score: PSM score obtained from Percolator
- q-value: q-value obtained from Percolator
- posterior_error_prob: posterior error probability obtained from Percolator

SD4: list of identified amino acid substitutions, tab-separated
- ProteinVariant: amino acid substitution description 
- Peptide: list of peptide sequences where the substitution has been identified, separated by semicolon
- PeptideType: type of the respected peptide
- PSMId: IDs of PSMs where this substitution has been identified
- Samples: IDs of respective samples

# References
\[1\] Spooner, W., McLaren, W., Slidel, T., Finch, D.K., Butler, R., Campbell, J., Eghobamien, L., Rider, D., Kiefer, C.M., Robinson, M.J., et al. (2018) Haplosaurus computes protein haplotypes for use in precision drug design. Nat. Commun., 9, 4128.
