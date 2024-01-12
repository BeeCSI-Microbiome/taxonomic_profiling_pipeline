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
ALL_FILES = [s+"_R1" for s in SAMPLE.index] + [s+"_R2" for s in SAMPLE.index]

# Create list of expected Kraken output file names
KRAKEN_REPORTS = expand('results/kraken/{sample}_report.txt', sample=SAMPLE.index)

# Krakefaction subworkflow variables
RAREFACTION_TABLES = expand("results/rarefied/{sample}.csv", sample=SAMPLE.index)

F_STRING = get_filter_string()


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#               Rules               #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# Define end goal output
rule all:
    input:
        KRAKEN_REPORTS, "results/multiqc_report.html", 'results/kronaplot.html',
        RAREFACTION_TABLES if config["perform_krakefaction"] else []
            
rule fastp:
    input:
        fwd=lambda wildcards: SAMPLE.loc[wildcards.sample, 'forward'],
        rev=lambda wildcards: SAMPLE.loc[wildcards.sample, 'reverse']
    output:
        fwd= retain(config["keep_fastp"], 'results/fastp/{sample}_fastp_R1.fastq.gz'),
        rev= retain(config["keep_fastp"], 'results/fastp/{sample}_fastp_R2.fastq.gz'),
        html= 'results/fastp/{sample}_fastp.html',
        json= 'results/fastp/{sample}_fastp.json'
    threads: 8
    conda:
        'envs/fastp.yaml'
    shell:
        'fastp -i {input.fwd} -I {input.rev} -o {output.fwd} -O {output.rev} --html {output.html} --json {output.json} --thread {threads}'

rule bowtie2:
    input:
        fwd= 'results/fastp/{sample}_fastp_R1.fastq.gz',
        rev= 'results/fastp/{sample}_fastp_R2.fastq.gz'
    output:
        fwd= retain(config["keep_bowtie2"], 'results/unmapped/{sample}_R1_unmapped.fastq.gz'),
        rev= retain(config["keep_bowtie2"], 'results/unmapped/{sample}_R2_unmapped.fastq.gz')
    threads: 8
    conda:
        'envs/bowtie2.yaml'
    log:
        'logs/bowtie2/{sample}.log'
    shell:
        '(bowtie2 -p {threads} -x phiX -1 {input.fwd} -2 {input.rev} --un-conc-gz results/unmapped/{wildcards.sample}_R%_unmapped.fastq.gz) > /dev/null 2> {log}'

rule fastqc:
    # wildcard 'readfile' is used because we must run fastqc on forward and reverse reads (R1, R2)
    input: 
        'results/unmapped/{readfile}_unmapped.fastq.gz'
    output:
        html= retain(config["keep_fastqc"], 'results/fastqc/{readfile}_unmapped_fastqc.html'),
        zipp= retain(config["keep_fastqc"], 'results/fastqc/{readfile}_unmapped_fastqc.zip')
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc {input} --outdir=results/fastqc'

rule multiqc:
    input: expand('results/fastqc/{readfile}_unmapped_fastqc.zip', readfile=ALL_FILES)
    output: 'results/multiqc_report.html'
    conda:
        'envs/multiqc.yaml'
    shell:
        "multiqc -o results -n multiqc_report.html {input}"

rule kraken2:
    input:
        fwd= 'results/unmapped/{sample}_R1_unmapped.fastq.gz',
        rev= 'results/unmapped/{sample}_R2_unmapped.fastq.gz'
    params:
        confidence = 0,
        base_qual = 0,
        db = config["db"]
    output:
        kraken_class = retain(config["keep_kraken_class"], 'results/kraken/{sample}_classification.txt'),
        kraken_report = 'results/kraken/{sample}_report.txt'
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

rule krona:
    input: expand('results/kraken/{sample}_report.txt', sample = SAMPLE.index)
    output:
        'results/kronaplot.html'
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
        'results/kraken/{sample}_classification.txt'
    output:
        retain(False, "results/filtered/{sample}_filtered.txt")
    params:
        filter_string=F_STRING
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

# custom translate script is used to produce the translated read files needed
# for krakefaction
rule translate:
    input:
        db_inspection = "results/db_inspection.txt",
        readfiles = expand("results/filtered/{sample}_filtered.txt", sample=SAMPLE.index)
    output:
        retain(False, expand("results/translated/{sample}_translated.txt", sample=SAMPLE.index))
    script:
        "scripts/kraken2-translate.py"

# perform rarefaction
rule krakefaction:
    input:
        trans = "results/translated/{sample}_translated.txt",
        untrans = "results/filtered/{sample}_filtered.txt"
    output:
        "results/rarefied/{sample}.csv"
    shell:
        "krakefaction -u {input.untrans} -t {input.trans} -o {output}"