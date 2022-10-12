configfile: "config.yaml"

rule all:
    input:
        in1="data/VariantDiscoverabilityList.tsv"

rule create_peptide_db:
    input:
        fasta=config['proteindb_fasta_stop'],
        subst_list="data/substitutions_list5.tsv"
    output:
        "data/PeptideList.tsv"
    params:
        max_cores=config['max_cores']
    threads: config['max_cores']
    shell:
        "python3 src/create_peptide_db.py -i {input.fasta} -m 2 -sl {input.subst_list} -t {params.max_cores} -o {output}"

rule db_aggregate_dupliate_peptides:
	input:
		"data/PeptideList.tsv"
	output:
		"data/PeptideListUniq.tsv"
	shell:
		"python3 src/check_duplicate_peptides.py -i {input} -o {output}"

rule db_get_multivar_peptides:
	input:
		"data/PeptideListUniq.tsv"
	output:
		"data/VariantPeptides.tsv"
	shell:
		"head -1 {input} > {output} ; grep \"variant\" {input} | grep -v \"reference\" >> {output}"

rule db_get_haplotype_freq:
	input:
		in1="data/VariantPeptides.tsv",
		in2="protein_haplotypes_1.csv"
	output:
		"data/VariantPeptidesFreq.tsv"
	shell:
		"python3 src/add_haplotype_frequency.py -i {input.in1} -hap {input.in2} -o {output}"

rule db_get_snp_discoverability:
    input:
        "data/VariantPeptidesFreq.tsv"
    output:
        "data/VariantDiscoverabilityList.tsv"
    shell:
        "python3 src/create_SNP_discoverability_list.py -i {input} -o {output}"

# ---------------- PSM post-processing ------------------

rule psm_annotate_variation:
    input:
        in1=config['merged_psm_list'],
        in2=config['proteindb_fasta'],
        in3="protein_haplotypes_1.csv"
    output:
        "data/PSM_reports_annotated.txt"
    params:
        max_cores=config['max_cores']
    threads: config['max_cores']
    shell:
        "python3 src/psm_check_haplotypes.py -i {input.in1} -f {input.in2} -hap {input.in3} -s all -t {params.max_cores} -o {output.out1}"

rule psm_threshold_FDR:
    input:
        "data/PSM_reports_annotated.txt"
    output:
        "data/PSM_reports_annotated_05FDR.txt"
    shell:
        "python3 src/psm_threshold.py -i {input} -thr_field q-value -thr_value 0.005 -o {output} "