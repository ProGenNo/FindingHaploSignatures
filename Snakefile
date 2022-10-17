configfile: "config.yaml"

rule all:
    input:
        db1="results/VariantDiscoverabilityList.tsv",
        db2="results/PeptideCoverageText.tsv",
        db3="results/PeptideCoverageCategoriesText.txt",
        psm1="results/PSM_protein_coverage_stats.tsv",
        psm2="results/PEP_violinplot.pdf",
        psm3="results/q-val_violinplot.pdf",
        psm4="results/ang_simil_violinplot.pdf",
        psm5="results/RT_diff_violinplot.pdf",
        psm6="results/PSMs_multivar_observed_predicted.txt",
        psm7="results/PSM_identified_variants.tsv",
        psm8="results/PSMs_gene_IDs.tsv"

rule create_peptide_db:
    input:
        fasta=config['proteindb_fasta_stop'],
        subst_list="data/substitutions_list5.tsv"
    output:
        "results/PeptideList.tsv"
    params:
        max_cores=config['max_cores']
    threads: config['max_cores']
    shell:
        "python3 src/create_peptide_db.py -i {input.fasta} -m 2 -sl {input.subst_list} -t {params.max_cores} -o {output}"

rule db_aggregate_dupliate_peptides:
	input:
		"results/PeptideList.tsv"
	output:
		all="results/PeptideListUniq.tsv",
		variant="results/VariantPeptides.tsv"
	shell:
		"python3 src/check_duplicate_peptides.py -i {input} -o {output.all} -var {output.variant}"

rule db_fill_missing_columns:
    input:
        uniq="results/PeptideListUniq.tsv",
        by_prot="results/PeptideList.tsv"
    output:
        "results/PeptideListComplete.tsv"
    params:
        max_cores=int(config['max_cores'] / 2)
    threads: config['max_cores'] / 2
    shell:
        "python3 src/db_fill_peptide_list_columns.py -uniq {input.uniq} -nonuniq {input.by_prot} -t {params.max_cores} -o {output}"

rule filter_fasta_stop:
    input:
        config['proteindb_fasta_stop']
    output:
        config['proteindb_fasta_stop_filtered']
    shell:
        "grep -A1 \"generic_ensref\" {input} > {output}; grep -A1 \"generic_enshap\" {input} >> {output}"

rule db_get_coverage_stats:
    input:
        pl="results/PeptideListComplete.tsv",
        fasta=config['proteindb_fasta_stop_filtered']
    output:
        "results/PeptideCoverageStats.tsv"
    params:
        max_cores=int(config['max_cores'] / 2)
    threads: config['max_cores'] / 2
    shell:
        "python3 src/get_peptide_stats_parallel.py -i {input.pl} -f {input.fasta}  -t {params.max_cores} -o {output} "    

rule db_get_coverage_categories:
    input:
        pl="results/PeptideListComplete.tsv",
        fasta=config['proteindb_fasta_stop_filtered']
    output:
        stats="results/PeptideCoverageCategories.tsv",
        text="results/PeptideCoverageCategoriesText.txt"
    params:
        max_cores=int(config['max_cores'] / 2)
    threads: config['max_cores'] / 2
    shell:
        "python3 src/get_peptide_stats_categories.py -i {input.pl} -f {input.fasta}  -t {params.max_cores} -o {output.stats} > {output.text} "  

rule db_store_coverage_stats_text:
    input:
        stats="results/PeptideCoverageStats.tsv",
        fasta=config['proteindb_fasta_stop_filtered']
    output:
        "results/PeptideCoverageText.tsv"
    shell:
        "python3 src/get_peptide_stats_precomputed.py -i {input.stats} -f {input.fasta} > {output}"

rule db_get_haplotype_freq:
    input:
        in1="results/VariantPeptides.tsv",
        in2="protein_haplotypes_1.csv"
    output:
        "results/VariantPeptidesFreq.tsv"
    params:
        max_cores=int(config['max_cores'] / 2)
    threads: config['max_cores'] / 2
    shell:
        "python3 src/add_haplotype_frequency.py -i {input.in1} -hap {input.in2} -t {params.max_cores} -o {output}"

rule db_get_snp_discoverability:
    input:
        "results/VariantPeptidesFreq.tsv"
    output:
        "results/VariantDiscoverabilityList.tsv"
    shell:
        "python3 src/create_SNP_discoverability_list.py -i {input} -o {output}"

# ---------------- PSM post-processing ------------------

rule psm_annotate_variation:
    input:
        in1=config['merged_psm_list'],
        in2=config['proteindb_fasta'],
        in3="protein_haplotypes_1.csv"
    output:
        "results/PSM_reports_annotated.txt"
    params:
        max_cores=config['max_cores']
    threads: config['max_cores']
    shell:
        "python3 src/psm_check_haplotypes.py -i {input.in1} -f {input.in2} -hap {input.in3} -s all -t {params.max_cores} -o {output}"

rule psm_threshold_FDR:
    input:
        "results/PSM_reports_annotated.txt"
    output:
        "results/PSM_reports_annotated_1FDR.txt"
    shell:
        "python3 src/psm_threshold.py -i {input} -thr_field q-value -thr_value 0.01 -o {output} "

rule filter_fasta:
    input:
        config['proteindb_fasta']
    output:
        config['proteindb_fasta_filtered']
    shell:
        "grep -A1 \"generic_ensref\" {input} > output; grep -A1 \"generic_enshap\" {input} >> {output}"

rule psm_protein_coverage_report:
    input:
        in1="results/PSM_reports_annotated_1FDR.txt",
        in2=config['proteindb_fasta_stop_filtered']
    output:
        "results/PSM_protein_coverage_report.tsv"
    shell:
        "python3 src/psm_check_coverage.py -i {input.in1} -f {input.in2} -o {output}"

rule psm_protein_coverage_stats:
    input:
        pep="results/PeptideCoverageStats.tsv",
        psm="results/PSM_protein_coverage_report.tsv",
        fasta=config['proteindb_fasta_stop_filtered']
    output:
        "results/PSM_protein_coverage_stats.tsv"
    shell:
        "python3 src/psm_get_coverage_stats.py -pep {input.pep} -psm {input.psm} -f {input.fasta} > {output} "

rule psm_assign_gene_id:
    input:
        psm="results/PSM_reports_annotated_1FDR.txt",
        fasta=config['proteindb_fasta_stop_filtered']
    output:
        "results/PSMs_gene_IDs.tsv"
    shell:
        "python src/psm_check_multigene.py -i {input.psm} -f {input.fasta} -o {output}"

rule psm_get_identified_variants:
    input:
        "results/PSM_reports_annotated_1FDR.txt"
    output:
        "results/PSM_identified_variants.tsv"
    shell:
        "python src/psm_identified_variants.py -i {input} -o {output}"

rule psm_get_observed_predicted_spectra:
    input:
        "results/PSM_reports_annotated_1FDR.txt"
    output:
        "results/PSMs_multivar_observed_predicted.txt"
    params:
        mzml_dir=config['mzml_dir'],
        predictions_dir=config['predictions_dir'],
        pep_type='multi-variant'
    shell:
        "python3 src/psm_get_observed_predicted_list.py -i {input} -m {params.mzml_dir} -p {params.predictions_dir} -type {params.pep_type} -o {output}"

rule psm_make_violinplots:
    input:
        "results/PSM_reports_annotated.txt"
    output:
        out1="results/PEP_violinplot.pdf",
        out2="results/q-val_violinplot.pdf",
        out3="results/ang_simil_violinplot.pdf",
        out4="results/RT_diff_violinplot.pdf"
    shell:
        "python3 src/psm_get_violinplots.py -i {input} -val_field posterior_error_prob -thr_field q-value -thr_value 0.01 -title \"PSM posterior error probability by peptide class\" -ylabel \"posterior error probability\" -o {output.out1} ; "
        "python3 src/psm_get_violinplots.py -i {input} -val_field q-value -thr_field q-value -thr_value 0.01 -title \"PSM q-value by peptide class\" -ylabel \"q-value\" -o {output.out2} ; "
        "python3 src/psm_get_violinplots.py -i {input} -val_field spectra_angular_similarity -thr_field q-value -thr_value 0.01 -title \"Spectrum similarity (predicted x observed) by peptide class\" -ylabel \"angular similarity\" -o {output.out3} ; "
        "python3 src/psm_get_violinplots.py -i {input} -val_field rt_Abs_error -thr_field q-value -thr_value 0.01 -title \"Retention time (RT) error by peptide class\" -ylabel \"RT error (sec.)\" -o {output.out4} ; "
