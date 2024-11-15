# make sure 1kg is using the same ref genome and dbsnp build annotation beforehand
# mapping file to update bim from grch37 to hrch38 _if needed_ - require specify ref genome beforehand?
# plink read/write w/o changes to correct ordering if needed
# annotate rsIDs w/ bcftools - do regardless
# have a few optional filters, like autosomal, biallelic etc.
# currently missing a sex check
# currently missing an INFO check/filter
# add rm-dup for samples
# add initial filter for v. extreme samples (+5 s.d.?)

PGEN_EXT = ["pvar", "psam", "pgen"]
SNPQC_OUT = ["afreq", "vmiss", "hardy", "log"]
SAMPLEQC_OUT = ["smiss", "het", "log"]
PRUNE_EXT = ["prune.in", "prune.out", "log"]
IBD_EXT = ["ibd", "log"]

configfile: "config.yaml"

rule get_snp_qc_stats:
    input:
        expand("raw_data/{input}.{ext}", input=config["input"]["plink_file"], ext=PGEN_EXT)
    params:
        infile=f"raw_data/{config['input']['plink_file']}",
        outfile="qc_info/01_snp_qc/snp_qc"
    output:
        expand("qc_info/01_snp_qc/snp_qc.{ext}", ext=SNPQC_OUT)
    threads:
        config["run"]["threads"]
    shell:
        """
        /home/matthew.smith/workdir2024/plink2 \
            --pfile {params.infile} \
            --freq \
            --hardy \
            --missing variant-only \
            --genotyping-rate \
            --out {params.outfile}
        """

# would be good if this updates based on snp QC somehow
# needs to handle INFO field and decide on multi-allelic or non-ATCG SNPs
rule filter_snps_by_qc:
    input:
        expand("raw_data/{input}.{ext}", input=config["input"]["plink_file"], ext=PGEN_EXT),
        expand("qc_info/01_snp_qc/snp_qc.{ext}", ext=SNPQC_OUT)
    params:
        infile=f"raw_data/{config['input']['plink_file']}",
        outfile="qc_info/02_post_snp_qc/qcd_data",
        maf=config['qc_thresholds']['maf'],
        mac=config['qc_thresholds']['mac'],
        hwe=config['qc_thresholds']['hwe'],
        vmiss=config['qc_thresholds']['vmiss'] 
    output:
        expand("qc_info/02_post_snp_qc/qcd_data.{ext}", ext=PGEN_EXT)
    threads:
        config["run"]["threads"]
    shell:
        """
        /home/matthew.smith/workdir2024/plink2 \
            --pfile {params.infile} \
            --autosome \
            --maf {params.maf} \
            --mac {params.mac} \
            --hwe {params.hwe} \
            --geno {params.vmiss} \
            --rm-dup force-first list \
            --make-pgen \
            --out {params.outfile}
        """

rule get_sample_qc_stats:
    input:
        expand("qc_info/02_post_snp_qc/qcd_data.{ext}", ext=PGEN_EXT)
    params:
        infile="qc_info/02_post_snp_qc/qcd_data",
        outfile="qc_info/03_het_smiss_sample_qc/sample_qc"
    output:
        expand("qc_info/03_het_smiss_sample_qc/sample_qc.{ext}", ext=SAMPLEQC_OUT)
    threads:
        config["run"]["threads"]
    shell:
        """
        /home/matthew.smith/workdir2024/plink2 \
            --pfile {params.infile} \
            --het \
            --missing sample-only \
            --out {params.outfile}
        """

# fix arg parsing in python script
rule find_heterozygosity_outliers:
    input:
        expand("qc_info/03_het_smiss_sample_qc/sample_qc.{ext}", ext=SAMPLEQC_OUT)
    params:
        infile="qc_info/02_post_snp_qc/qcd_data.het",
        het=config['qc_thresholds']['het']
    output:
        "qc_info/03_het_smiss_sample_qc/het_outliers.txt"
    threads:
        config["run"]["threads"]
    shell:
        """
        python3 scripts/filter_ids_by_het.py --infile {params.infile} --n_sds {params.het} --outfile {output}
        """

rule filter_samples_by_qc:
    input:
        "qc_info/03_het_smiss_sample_qc/het_outliers.txt",
        expand("qc_info/02_post_snp_qc/qcd_data.{ext}", ext=PGEN_EXT)
    params:
        infile="qc_info/02_post_snp_qc/qcd_data",
        filterfile="qc_info/03_het_smiss_sample_qc/het_outliers.txt",
        outfile="qc_info/04_post_het_smiss_sample_qc/qcd_data",
        smiss=config['qc_thresholds']['smiss']
    output:
        expand("qc_info/04_post_het_smiss_sample_qc/qcd_data.{ext}", ext=PGEN_EXT)
    threads:
        config["run"]["threads"]
    shell:
        """
        /home/matthew.smith/workdir2024/plink2 \
            --pfile {params.infile} \
            --remove {params.filterfile} \
            --mind {params.smiss} \
            --remove-nosex \
            --make-pgen \
            --out {params.outfile}
        """

rule prune_for_qc:
    input:
        expand("qc_info/04_post_het_smiss_sample_qc/qcd_data.{ext}", ext=PGEN_EXT)
    params:
        infile="qc_info/04_post_het_smiss_sample_qc/qcd_data",
        pihat=config['qc_thresholds']['pihat'],
        window=config['qc_thresholds']['qc_prune_window'],
        step=config['qc_thresholds']['qc_prune_step'],
        r2=config['qc_thresholds']['qc_prune_r2'],
        outfile="qc_info/05_relatedness_ancestry_checks/qc_prune"
    output:
        expand("qc_info/05_relatedness_ancestry_checks/qc_prune.{ext}", ext=PRUNE_EXT)
    threads:
        config["run"]["threads"]
    shell:
        """
        /home/matthew.smith/workdir2024/plink2 \
            --pfile {params.infile} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --out {params.outfile}
        """

# decide on plink1 with --genome or plink2 if king table
# rule check_relatedness:
#     input:
#         expand("qc_info/05_relatedness_ancestry_checks/qc_prune.{ext}", ext=PRUNE_EXT),
#         expand("qc_info/04_post_het_smiss_sample_qc/qcd_data.{ext}", ext=PGEN_EXT)
#     params:
#         infile="qc_info/04_post_het_smiss_sample_qc/qcd_data",
#         filterfile="qc_info/05_relatedness_ancestry_checks/qc_prune.prune.in",
#         maf=config['qc_thresholds']['qc_prune_maf'],
#         outfile="qc_info/05_relatedness_ancestry_checks/relatedness_check"
#     output:
#         expand("qc_info/05_relatedness_ancestry_checks/relatedness_check.{ext}", ext=IBD_EXT)
#     shell:
#         """
#         /home/matthew.smith/workdir2024/plink2 \
#             --pfile {params.infile} \
#             --extract {params.filterfile} \
#             --maf 0.05 \
#             --genome \
#             --out 


# rule find_related_individuals:

# rule filter_samples_by_relatedness:

# rule calculate_principal_components:

# optional restrict to population by merging w/ 1kg 

# rule regenerating_snp_summaries

# rule regenerate_sample_summaries

# rule generate_summary_report

# rule all:
# input:
# "my_qcd_files.bed/bim/fam or /pgen/pvar/psam

# randomly subset samples into train and test

# generate report on covariates etc. by train/test split

# export as additive SNPs with reference allele specified for train/test

# shuffle .raw file with shuf

# read raw with numpy and convert to hdf5 with dask chunks

# run final check that hdf5 files read correctly and have no overlap in samples but complete overlap in SNPs for train/test