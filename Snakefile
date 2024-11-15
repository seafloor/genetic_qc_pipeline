# make sure 1kg is using the same ref genome and dbsnp build annotation beforehand
# mapping file to update bim from grch37 to hrch38 _if needed_ - require specify ref genome beforehand?
# plink read/write w/o changes to correct ordering if needed
# annotate rsIDs w/ bcftools - do regardless
# have a few optional filters, like autosomal, biallelic etc.
# currently missing a sex check
# currently missing an INFO check/filter
# add rm-dup for samples

PGEN_EXT = ["pvar", "psam", "pgen"]
SNPQC_OUT = ["afreq", "vmiss", "hardy", "log"]
SAMPLEQC_OUT = ["smiss", "het", "log"]

configfile: "config.yaml"

rule get_snp_qc_stats:
    input:
        expand("raw_data/{input}.{ext}", input=config["input"]["plink_file"], ext=PGEN_EXT)
    params:
        infile=f"raw_data/{config['input']['plink_file']}",
        outfile="qc_info/01_snp_qc/snp_qc"
    output:
        expand("qc_info/01_snp_qc/snp_qc.{ext}", ext=SNPQC_OUT)
    threads: config["run"]["threads"]
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
    threads: config["run"]["threads"]
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
        outfile="qc_info/03_sample_qc/sample_qc"
    output:
        expand("qc_info/03_sample_qc/sample_qc.{ext}", ext=SAMPLEQC_OUT)
    threads: config["run"]["threads"]
    shell:
        """
        /home/matthew.smith/workdir2024/plink2 \
            --pfile {params.infile} \
            --het \
            --missing sample-only \
            --out {params.outfile}
        """

rule filter_samples_by_qc:
    input:
        expand("qc_info/03_sample_qc/sample_qc.{ext}", ext=SAMPLEQC_OUT)
    params:
        infile="qc_info/02_post_snp_qc/qcd_data.het",
        outfile=outfile="qc_info/03_post_sample_qc/qcd_data",
        het="",
        smiss=config['qc_thresholds']['smiss']
    output:
        expand("qc_info/03_post_sample_qc/qcd_data.{ext}", ext=PGEN_EXT)
    threads: config["run"]["threads"]
    shell:
        """
        python3 scripts/filter_ids_by_het.py --infile --n_sds --outfile
        # awk command to filter IDs from het file to ID file > het_outlier_ids.txt
        /home/matthew.smith/workdir2024/plink2 \
            --pfile {params.infile} \
            --remove het_outlier_ids.txt \
            --mind {params.smiss} \
            --remove-nosex \
            --make-pgen \
            --out {params.outfile}

# rule filter_samples_by_qc
# include remove nosex, rm-dup

# optional restrict to population by merging w/ 1kg 

# rule regenerating_snp_summaries

# rule regenerate_sample_summaries

# rule generate_summary_report

# rule all:
# input:
# "my_qcd_files.bed/bim/fam or /pgen/pvar/psam
