import pandas as pd

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#               Setup               #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# Point to config file in which sample sheet, output path and kraken database are specified
configfile: 'config/config_kaiju.yaml'

# Get functions
include: "rules/common.smk"

# Get sample prefixes
SAMPLE = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
ALL_FILES = [s+"_R1" for s in SAMPLE.index] + [s+"_R2" for s in SAMPLE.index]

TAXA_RANKS = ["phylum","class","order","family","genus","species"]


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#               Rules               #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# Define end goal output
rule all:
    input:
        "results/multiqc_report.html",
        expand('results/kaiju_summaries/kaiju_{rank}_summary.tsv', rank=TAXA_RANKS),
        expand('results/krona/{sample}_kronaplot.html', sample=SAMPLE.index)

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
    log: 'logs/fastp/fastp_{sample}.out'
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

rule kaiju:
    input:
        fwd= 'results/unmapped/{sample}_R1_unmapped.fastq.gz',
        rev= 'results/unmapped/{sample}_R2_unmapped.fastq.gz'
    params:
        DB = config["db"],
        NODES_DMP = config["NODES_DMP"]
    output:
        kaiju_out =  'results/kaiju/{sample}.kaiju.out',
    threads: 16
    resources:
        mem_mb=256000,
        disk_mb=256000,
        res_tmpfs=256000
    conda:
        'envs/kaiju.yaml'
    shell:
        "kaiju "
        "-t {params.NODES_DMP} "
        "-f {params.DB} "
        "-z {threads} "
        "-o {output.kaiju_out} "
        "-i {input.fwd} "
        "-j {input.rev}"

rule kaiju2krona:
    input: 'results/kaiju/{sample}.kaiju.out'
    output: 'results/kaiju2krona/{sample}.kaiju.out.krona'
    params:
        NODES_DMP = config["NODES_DMP"],
        NAMES_DMP = config["NAMES_DMP"]
    conda:
        'envs/kaiju.yaml'
    shell:
        "kaiju2krona "
        "-t {params.NODES_DMP} "
        "-n {params.NAMES_DMP} "
        "-i {input} "
        "-o {output} "
        "-u"

rule krona:
    input: 'results/kaiju2krona/{sample}.kaiju.out.krona'
    output: 'results/krona/{sample}_kronaplot.html'
    conda:
        'envs/krona.yaml'
    shell:
        "ktImportText -o {output} {input}"

rule kaiju2table_species:
    input: expand('results/kaiju/{sample}.kaiju.out', sample = SAMPLE.index)
    output:
        SPECIES='results/kaiju_summaries/kaiju_species_summary.tsv'
    params:
        NODES_DMP = config["NODES_DMP"],
        NAMES_DMP = config["NAMES_DMP"]
    conda:
        'envs/kaiju.yaml'
    resources:
        mem_mb=256000,
        disk_mb=256000,
        res_tmpfs=256000
    shell:
        "kaiju2table -t {params.NODES_DMP} -n {params.NAMES_DMP} -p -r species -o {output.SPECIES} {input}"

rule kaiju2table_genus:
    input: expand('results/kaiju/{sample}.kaiju.out', sample = SAMPLE.index)
    output: 
        GENUS='results/kaiju_summaries/kaiju_genus_summary.tsv',
    params:
        NODES_DMP = config["NODES_DMP"],
        NAMES_DMP = config["NAMES_DMP"]
    conda:
        'envs/kaiju.yaml'
    resources:
        mem_mb=256000,
        disk_mb=256000,
        res_tmpfs=256000
    shell:
        "kaiju2table -t {params.NODES_DMP} -n {params.NAMES_DMP} -p -r genus -o {output.GENUS} {input}"

rule kaiju2table_family:
    input: expand('results/kaiju/{sample}.kaiju.out', sample = SAMPLE.index)
    output: 
        FAMILY='results/kaiju_summaries/kaiju_family_summary.tsv',
    params:
        NODES_DMP = config["NODES_DMP"],
        NAMES_DMP = config["NAMES_DMP"]
    conda:
        'envs/kaiju.yaml'
    resources:
        mem_mb=256000,
        disk_mb=256000,
        res_tmpfs=256000
    shell:
        "kaiju2table -t {params.NODES_DMP} -n {params.NAMES_DMP} -p -r family -o {output.FAMILY} {input}"

rule kaiju2table_order:
    input: expand('results/kaiju/{sample}.kaiju.out', sample = SAMPLE.index)
    output: 
        ORDER='results/kaiju_summaries/kaiju_order_summary.tsv',
    params:
        NODES_DMP = config["NODES_DMP"],
        NAMES_DMP = config["NAMES_DMP"]
    conda:
        'envs/kaiju.yaml'
    resources:
        mem_mb=256000,
        disk_mb=256000,
        res_tmpfs=256000
    shell:
        "kaiju2table -t {params.NODES_DMP} -n {params.NAMES_DMP} -p -r order -o {output.ORDER} {input}"
        
rule kaiju2table_class:
    input: expand('results/kaiju/{sample}.kaiju.out', sample = SAMPLE.index)
    output: 
        CLASS='results/kaiju_summaries/kaiju_class_summary.tsv',
    params:
        NODES_DMP = config["NODES_DMP"],
        NAMES_DMP = config["NAMES_DMP"]
    conda:
        'envs/kaiju.yaml'
    resources:
        mem_mb=256000,
        disk_mb=256000,
        res_tmpfs=256000
    shell:
        "kaiju2table -t {params.NODES_DMP} -n {params.NAMES_DMP} -p -r class -o {output.CLASS} {input}"

rule kaiju2table_phylum:
    input: expand('results/kaiju/{sample}.kaiju.out', sample = SAMPLE.index)
    output: 
        PHYLUM='results/kaiju_summaries/kaiju_phylum_summary.tsv',
    params:
        NODES_DMP = config["NODES_DMP"],
        NAMES_DMP = config["NAMES_DMP"]
    conda:
        'envs/kaiju.yaml'
    resources:
        mem_mb=256000,
        disk_mb=256000,
        res_tmpfs=256000
    shell:
        "kaiju2table -t {params.NODES_DMP} -n {params.NAMES_DMP} -p -r phylum -o {output.PHYLUM} {input}"
