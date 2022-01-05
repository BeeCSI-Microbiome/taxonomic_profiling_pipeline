import pandas as pd

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#               Setup               #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# Point to config file in which sample sheet, output path and kraken database are specified
configfile: 'config.yaml'

# Get output path
OUTDIR = config["outdir"]
# Get sample prefixes
SAMPLE = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
ALL_FILES = [s+"_R1" for s in SAMPLE.index] + [s+"_R2" for s in SAMPLE.index]

# Get filter target taxa
FILTER_TARGETS = config['filter_targets']

# Create list of expected Kraken output file names
KRAKEN_REPORTS = expand(OUTDIR + '/kraken/{sample}_report.txt', sample=SAMPLE.index)


# Krakefaction subworkflow variables
RAREFACTION_TABLES = expand(OUTDIR + "/rarefied/{sample}.csv", sample=SAMPLE.index)
# create the filter string based on config file
F_STRING="unclassified\|cellular organism"
if FILTER_TARGETS:
	F_STRING = F_STRING + "\|" + "\|".join(FILTER_TARGETS)
print("Filtering with the following filter string:\n\t'{}'".format(F_STRING))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#             Functions             #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

def retain(flag, path):
    """Returns given path if flag is true, else returns temp(path)
       which causes the output to be deleted after it is no longer needed"""
    if flag: return path
    else: return temp(path)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#               Rules               #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# Define end goal output
rule all:
    input:
        KRAKEN_REPORTS, OUTDIR+"/multiqc_report.html", OUTDIR + '/kronaplot.html',
        RAREFACTION_TABLES if config["perform_krakefaction"] else []
            

rule fastp:
    input:
        fwd=lambda wildcards: SAMPLE.loc[wildcards.sample, 'forward'],
        rev=lambda wildcards: SAMPLE.loc[wildcards.sample, 'reverse']
    output:
        fwd= retain(config["keep_fastp"], OUTDIR + '/fastp/{sample}_fastp_R1.fastq.gz'),
        rev= retain(config["keep_fastp"], OUTDIR + '/fastp/{sample}_fastp_R2.fastq.gz'),
        html= OUTDIR + '/fastp/{sample}_fastp.html',
        json= OUTDIR + '/fastp/{sample}_fastp.json'
    threads: 8
    conda:
        'envs/fastp.yaml'
    shell:
        'fastp -i {input.fwd} -I {input.rev} -o {output.fwd} -O {output.rev} --html {output.html} --json {output.json} --thread {threads}'

rule bowtie2:
    input:
        fwd= OUTDIR + '/fastp/{sample}_fastp_R1.fastq.gz',
        rev= OUTDIR + '/fastp/{sample}_fastp_R2.fastq.gz'
    output:
        fwd= retain(config["keep_bowtie2"], OUTDIR + '/unmapped/{sample}_R1_unmapped.fastq.gz'),
        rev= retain(config["keep_bowtie2"], OUTDIR + '/unmapped/{sample}_R2_unmapped.fastq.gz')
    params:
        out=OUTDIR
    threads: 8
    conda:
        'envs/bowtie2.yaml'
    log:
        'logs/bowtie2/{sample}.log'
    shell:
        '(bowtie2 -p {threads} -x phiX -1 {input.fwd} -2 {input.rev} --un-conc-gz {params.out}/unmapped/{wildcards.sample}_R%_unmapped.fastq.gz) 2> {log}'


rule fastqc:
    # wildcard 'readfile' is used because we must now run fastqc on forward and reverse reads
    input: 
        OUTDIR + '/unmapped/{readfile}_unmapped.fastq.gz'
    output:
        html= retain(config["keep_fastqc"], OUTDIR + '/fastqc/{readfile}_unmapped_fastqc.html'),
        zipp= retain(config["keep_fastqc"], OUTDIR + '/fastqc/{readfile}_unmapped_fastqc.zip')
    params:
        outd=OUTDIR
    conda:
        'envs/fastqc.yaml'
    shell:
        'fastqc {input} --outdir={params.outd}/fastqc'

rule multiqc:
    input: expand(OUTDIR + '/fastqc/{readfile}_unmapped_fastqc.zip', readfile=ALL_FILES)
    output: OUTDIR + '/multiqc_report.html'
    conda:
        'envs/multiqc.yaml'
    params:
        outd=OUTDIR
    shell:
        "multiqc -o {params.outd} -n multiqc_report.html {input}"

rule kraken2:
    input:
        fwd= OUTDIR + '/unmapped/{sample}_R1_unmapped.fastq.gz',
        rev= OUTDIR + '/unmapped/{sample}_R2_unmapped.fastq.gz'
    params:
        confidence = 0,
        base_qual = 0,
        db = config["db"]
    output:
        kraken_class = retain(config["keep_kraken_class"], OUTDIR + '/kraken/{sample}_classification.txt'),
        kraken_report = OUTDIR + '/kraken/{sample}_report.txt'
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
    input: expand(OUTDIR + '/kraken/{sample}_report.txt', sample = SAMPLE.index)
    output:
        OUTDIR + '/kronaplot.html'
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
        OUTDIR + '/kraken/{sample}_classification.txt'
    output:
        retain(False, OUTDIR + "/filtered/{sample}_filtered.txt")
    shell:
        "cat {input} | grep -v \"" + F_STRING + "\" > {output}"

# Get the db_inspection file needed for the translate step
rule get_db_inspection:
    input:
        DB = config["db"]
    output:
        OUTDIR + "/db_inspection.txt"
    conda:
        'envs/kraken2.yaml'
    shell:
        "kraken2-inspect --db {input} > {output}"

# custom translate script is used to produce the translated read files needed
# for krakefaction
rule translate:
    input:
        db_inspection = OUTDIR + "/db_inspection.txt",
        readfiles = expand(OUTDIR + "/filtered/{sample}_filtered.txt", sample=SAMPLE.index)
    output:
        retain(False, expand(OUTDIR + "/translated/{sample}_translated.txt", sample=SAMPLE.index))
    script:
        "scripts/kraken2-translate.py"

# perform rarefaction
rule krakefaction:
    input:
        trans = OUTDIR + "/translated/{sample}_translated.txt",
        untrans = OUTDIR + "/filtered/{sample}_filtered.txt"
    output:
        OUTDIR + "/rarefied/{sample}.csv"
    shell:
        "krakefaction -u {input.untrans} -t {input.trans} -o {output}"