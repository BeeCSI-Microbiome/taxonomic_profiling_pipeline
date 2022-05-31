import pandas as pd

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#               Setup               #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# Point to config file in which sample sheet, output path and kraken database are specified
configfile: 'config.yaml'

# Get functions
include: "rules/common.smk"

# Get sample prefixes
SAMPLE = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
# ALL_FILES = [s+"_R1" for s in SAMPLE.index] + [s+"_R2" for s in SAMPLE.index]
DATASET_NAMES = set(SAMPLE.dataset)

def get_kraken_reports():
    l = []
    for d in DATASET_NAMES:
        l.extend(expand('results/{ds}/kraken/{sample}_report.txt', ds=d, sample = SAMPLE.index[SAMPLE['dataset'] == d].tolist()))
    return l

# Krakefaction subworkflow variables
def get_krakefaction_tables():
    if config["perform_krakefaction"]:
        l = []
        for d in DATASET_NAMES:
            l.extend(expand('results/{ds}/rarefied/{sample}.csv', ds=d, sample = SAMPLE.index[SAMPLE['dataset'] == d].tolist()))
        return l
    else:
        return []

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#               Rules               #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# Define end goal output
rule all:
    input:
        get_kraken_reports(),
        expand("results/{dataset}/multiqc_report.html", dataset=DATASET_NAMES), 
        expand('results/{dataset}/kronaplot.html', dataset=DATASET_NAMES),
        get_krakefaction_tables()
            
rule fastp:
    input:
        fwd=lambda wildcards: SAMPLE.loc[wildcards.sample, 'forward'],
        rev=lambda wildcards: SAMPLE.loc[wildcards.sample, 'reverse']
    output:
        fwd= retain(config["keep_fastp"], 'results/{dataset}/fastp/{sample}_fastp_R1.fastq.gz'),
        rev= retain(config["keep_fastp"], 'results/{dataset}/fastp/{sample}_fastp_R2.fastq.gz'),
        html= 'results/{dataset}/fastp/{sample}_fastp.html',
        json= 'results/{dataset}/fastp/{sample}_fastp.json'
    threads: 16
    conda:
        'envs/fastp.yaml'
    shell:
        'fastp -i {input.fwd} -I {input.rev} -o {output.fwd} -O {output.rev} --html {output.html} --json {output.json} --thread {threads}'

rule bowtie2:
    input:
        fwd= 'results/{dataset}/fastp/{sample}_fastp_R1.fastq.gz',
        rev= 'results/{dataset}/fastp/{sample}_fastp_R2.fastq.gz'
    output:
        fwd= retain(config["keep_bowtie2"], 'results/{dataset}/unmapped/{sample}_R1_unmapped.fastq.gz'),
        rev= retain(config["keep_bowtie2"], 'results/{dataset}/unmapped/{sample}_R2_unmapped.fastq.gz')
    threads: 16
    conda:
        'envs/bowtie2.yaml'
    log:
        'logs/bowtie2/{dataset}/{sample}.log'
    shell:
        '(bowtie2 -p {threads} -x phiX -1 {input.fwd} -2 {input.rev} --un-conc-gz results/{wildcards.dataset}/unmapped/{wildcards.sample}_R%_unmapped.fastq.gz) > /dev/null 2> {log}'

def get_bowtie_checkpoint_input(wildcards):
    return expand('results/{{dataset}}/unmapped/{sample}_R{n}_unmapped.fastq.gz', sample = SAMPLE.index[SAMPLE['dataset'] == wildcards.dataset].tolist(), n = [1,2])

# Checkpoint to ensure bowtie is completed for a dataset before running subsequent rules, thus ensuring temp bowtie input (fastp output)
# is deleted, saving space for subsequent jobs
rule bowtie_checkpoint:
    input:
        get_bowtie_checkpoint_input
    output:
        'results/{dataset}/bowtie_complete'
    shell:
        'touch {output}'

rule fastqc:
    # wildcard 'readfile' is used because we must run fastqc on forward and reverse reads (R1, R2)
    input:
        'results/{dataset}/bowtie_complete',
        'results/{dataset}/unmapped/{readfile}_unmapped.fastq.gz'
    output:
        html= retain(config["keep_fastqc"], 'results/{dataset}/fastqc/{readfile}_unmapped_fastqc.html'),
        zipp= retain(config["keep_fastqc"], 'results/{dataset}/fastqc/{readfile}_unmapped_fastqc.zip')
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc {input} --outdir=results/{wildcards.dataset}/fastqc'

def get_multiqc_input(wildcards):
    return expand('results/{{dataset}}/fastqc/{sample}_R{n}_unmapped_fastqc.zip', sample = SAMPLE.index[SAMPLE['dataset'] == wildcards.dataset].tolist(), n = [1,2])

rule multiqc:
    input: 
        get_multiqc_input
    output: 'results/{dataset}/multiqc_report.html'
    conda:
        'envs/multiqc.yaml'
    shell:
        "multiqc -o results/{wildcards.dataset} -n multiqc_report.html {input}"

rule kraken2:
    input:
        'results/{dataset}/bowtie_complete',
        'results/{dataset}/multiqc_report.html',
        fwd= 'results/{dataset}/unmapped/{sample}_R1_unmapped.fastq.gz',
        rev= 'results/{dataset}/unmapped/{sample}_R2_unmapped.fastq.gz'
    params:
        confidence = 0,
        base_qual = 0,
        db = config["db"]
    output:
        kraken_class = retain(config["keep_kraken_class"], 'results/{dataset}/kraken/{sample}_classification.txt'),
        kraken_report = 'results/{dataset}/kraken/{sample}_report.txt'
    threads: 16
    conda:
        'envs/kraken2.yaml'
    shell:
        "kraken2 "
        "--db {params.db} "
        "--threads {threads} "
        "--output {output.kraken_class} "
        "--confidence {params.confidence} "
        "--minimum-base-quality {params.base_qual} "
        "--report {output.kraken_report} "
        "--paired "
        "--use-names "
        "--gzip-compressed "
        "{input.fwd} {input.rev}"

def get_kraken_checkpoint_input(wildcards):
    return expand('results/{{dataset}}/kraken/{sample}_report.txt', sample = SAMPLE.index[SAMPLE['dataset'] == wildcards.dataset].tolist())

# Checkpoint to ensure kraken is completed for a dataset before running subsequent rules, thus ensuring temp bowtie output
# is deleted, saving space for subsequent jobs
rule kraken_checkpoint:
    input:
        get_kraken_checkpoint_input
    output:
        'results/{dataset}/kraken_complete'
    shell:
        'touch {output}'

def get_krona_input(wildcards):
    # get filepaths expanded to contain only samples from the current dataset
    return expand('results/{{dataset}}/kraken/{sample}_report.txt', sample = SAMPLE.index[SAMPLE['dataset'] == wildcards.dataset].tolist())

rule krona:
    input:
        get_krona_input
    output:
        'results/{dataset}/kronaplot.html'
    conda:
        'envs/krona.yaml'
    shell:
        """ktUpdateTaxonomy.sh
        ktImportTaxonomy -m 3 -t 5 {input} -o {output}"""

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#           Krakefaction            #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# filter out unclassifed reads, "cellular organism" reads, and all others specified
# in the config file
rule filter:
    input:
        'results/{dataset}/kraken_complete',
        'results/{dataset}/kraken/{sample}_classification.txt'
    output:
        temp("results/{dataset}/filtered/{sample}_filtered.txt")
    params:
        filter_string=get_filter_string()
    shell:
        "cat {input} | grep -v \"{params.filter_string}\" > {output}"

# Get the db_inspection file needed for the translate step
rule get_db_inspection:
    input:
        DB = config["db"]
    output:
        "results/db_inspection.txt"
    conda:
        'envs/kraken2.yaml'
    shell:
        "kraken2-inspect --db {input} > {output}"

def get_krakefaction_translate_input():
    l = []
    for d in DATASET_NAMES:
        l.extend(expand('results/{ds}/filtered/{sample}_filtered.txt', ds=d, sample = SAMPLE.index[SAMPLE['dataset'] == d].tolist()))
    return l

def get_krakefaction_translate_output():
    l = []
    for d in DATASET_NAMES:
        l.extend(expand("results/{ds}/translated/{sample}_translated.txt", ds=d, sample = SAMPLE.index[SAMPLE['dataset'] == d].tolist()))
    return l

# custom translate script is used to produce the translated read files needed
# for krakefaction
rule translate:
    input:
        # db_inspection must be first argument to script
        db_inspection = "results/db_inspection.txt",
        input_files = get_krakefaction_translate_input()
    output:
        temp(get_krakefaction_translate_output())
    script:
        "scripts/kraken2-translate.py"

# perform rarefaction
rule krakefaction:
    input:
        trans = "results/{dataset}/translated/{sample}_translated.txt",
        untrans = "results/{dataset}/filtered/{sample}_filtered.txt"
    output:
        "results/{dataset}/rarefied/{sample}.csv"
    shell:
        "krakefaction -u {input.untrans} -t {input.trans} -o {output}"