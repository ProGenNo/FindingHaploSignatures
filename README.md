# Identifying Protein Haplotypes by Mass Spectrometry
Code related to the "Identifying Protein Haplotypes by Mass Spectrometry" publication

## Requirements and Usage
Required software is Snakemake and Conda, rest of libraries are included in the provided Conda environment, created automatically by Snakemake.

Steps for reproducing results:
- Download supplementary data from (LINK)
- Clone this repository
- Provide path to the downloaded files following instructions `config.yaml`
- Run the pipeline using `snakemake -c<# cores> -p --use-conda`

## Additional files provided
Two additional files were derived from data made public by Spooner et al. \[1\]

`protein_haplotypes_1.csv`: List of protein haplotypes created using alleles with minor allele frequency (MAF) > 0.01
`substitutions_list.tsv`: List of substitutions found in the provided haplotypes, along with other pre-computed values 

# References
\[1\] Spooner, W., McLaren, W., Slidel, T., Finch, D.K., Butler, R., Campbell, J., Eghobamien, L., Rider, D., Kiefer, C.M., Robinson, M.J., et al. (2018) Haplosaurus computes protein haplotypes for use in precision drug design. Nat. Commun., 9, 4128.
