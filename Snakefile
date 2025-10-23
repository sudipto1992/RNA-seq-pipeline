configfile: "config.yaml"

import glob
import os


MYGTF = config["gtf_file"]
SAMPLES = config["samples"]
print ("sample_loaded:", SAMPLES)

os.makedirs("log", exist_ok=True)
os.makedirs("genome_build", exist_ok=True)

# QC steps
target = []
if config["downloadSRR"]:
    target+=expand("fastq_input/{sample}_1.fastq",sample=SAMPLES),
    target+=expand("fastq_input/{sample}_2.fastq",sample=SAMPLES)
if config["iscompressed"]:
    target+=expand("fastq_input/{sample}_1.fastq",sample=SAMPLES),
    target+=expand("fastq_input/{sample}_2.fastq",sample=SAMPLES)
if config['runtrimming']:
    target+=expand("trimmed/{sample}_1_paired.fastq",sample=SAMPLES),
    target+=expand("trimmed/{sample}_2_paired.fastq",sample=SAMPLES)
if config["gen_fastqc"]:
    target+=expand("fastqc/{sample}_1_paired_fastqc.html",sample=SAMPLES),
    target+=expand("fastqc/{sample}_2_paired_fastqc.html",sample=SAMPLES)
if config["gen_multiqc"]:
    target.append("multiqc/multiqc_report.html")
if config["genome_anno"]:
    target.append("annotation/organism.ss"),
    target.append("annotation/organism.exon")
if config["genome_build"]:
    target.append(f"{config['genome_index_dir']}/genome.1.ht2")
if config["use_hisat"]:
    target+= expand("final/aligned_{sample}.bam", sample=SAMPLES)
if config["use_featurecount"]:
    target+= expand("final/{sample}_count.txt",sample=SAMPLES)
print (target)

# What is the final target?
rule all:
    input:
        target

# local rules for steps
rule downloadSRR:
    output:
        R1="fastq_input/{sample}_1.fastq",
        R2="fastq_input/{sample}_2.fastq"
    log:
        "log/datacol_{sample}.log"
    conda:
        "envs/QC.yaml"
    shell:
        """
        set -x
        mkdir -p sra_raw fastq_input
        prefetch {wildcards.sample} --output-directory sra_raw
        fasterq-dump sra_raw/{wildcards.sample}/{wildcards.sample}.sra --threads 8 --outdir fastq_input
        rm -r sra_raw/{wildcards.sample}
        """
rule is_compressed:
    input:
        P1="fastq_input/{sample}_1.fastq.gz",
        P2="fastq_input/{sample}_2.fastq.gz"
    output:
        R1="fastq_input/{sample}_1.fastq",
        R2="fastq_input/{sample}_2.fastq"
    log:
        "log/unzip_{sample}.log"
    shell:
        """
        set -x
        gunzip -c {input.P1} > {output.R1}
        gunzip -c {input.P2} > {output.R2}
        """

rule run_trimmomatic:
    input:
        R1="fastq_input/{sample}_1.fastq",
        R2="fastq_input/{sample}_2.fastq"
    output:
        R1_paired="trimmed/{sample}_1_paired.fastq",
        R1_unpaired="trimmed/{sample}_1_unpaired.fastq",
        R2_paired="trimmed/{sample}_2_paired.fastq",
        R2_unpaired="trimmed/{sample}_2_unpaired.fastq"
    params:
        adapters="TruSeq3-PE.fa",
        minlen=config["trimmotatic_minlen"]
    log:
        "log/trimmomatic_{sample}.log"
    conda:
        "envs/QC.yaml"
    shell:
        """
        mkdir -p trimmed
        trimmomatic PE -threads 4 \
            {input.R1} {input.R2} \
            {output.R1_paired} {output.R1_unpaired} \
            {output.R2_paired} {output.R2_unpaired} \
            ILLUMINACLIP:{params.adapters}:2:30:10\
            MINLEN:{params.minlen} \
            > {log} 2>&1
        rm -r {input.R1}
        rm -r {input.R2}
        """

    
rule gen_fastqc:
    input:
        P1="trimmed/{sample}_1_paired.fastq",
        P2="trimmed/{sample}_2_paired.fastq"
    output:
        R1="fastqc/{sample}_1_paired_fastqc.html",
        R2="fastqc/{sample}_2_paired_fastqc.html"
    log:
        "log/qc_{sample}.log"
    conda:
        "envs/QC.yaml"
    shell:
        """
        set -x
        mkdir -p fastqc
        fastqc {input.P1} {input.P2} -o fastqc > {log} 2>&1
        """

rule gen_multiqc:
    input:
        expand("fastqc/{sample}_1_paired_fastqc.html",sample=SAMPLES),
        expand("fastqc/{sample}_2_paired_fastqc.html",sample=SAMPLES)
    output:
        "multiqc/multiqc_report.html"
    log:
        "log/qc.log"
    conda:
        "envs/QC.yaml"
    shell:
        """
        set -x
        mkdir -p multiqc
        multiqc fastqc -o multiqc > {log} 2>&1
        """

rule genome_anno:
    input:
        config["gtf_file"]
    output:
        R1="annotation/organism.ss",
        R2="annotation/organism.exon"
    log:
        "log/genome_build.log"
    conda:
        "envs/rnaseq.yaml"
    shell:
        """
        set -x
        hisat2_extract_splice_sites.py {input} > {output.R1} > {log} 2>&1
        hisat2_extract_exons.py {input} > {output.R2} > {log} 2>&1
        """

rule genome_build:
    input:
        fasta = config["genome"],
        ss = "annotation/organism.ss",
        exon = "annotation/organism.exon"
    output:
        idx = config["genome_index_dir"]+"/genome.1.ht2"
    log:
        "log/genome_build.log"
    conda:
        "envs/rnaseq.yaml"
    shell:
        """
        set -x
        hisat2-build -p 8 --ss {input.ss} --exon {input.exon} {input.fasta} {config[genome_index_dir]}/genome > {log} 2>&1
        """

rule use_hisat:
    input:
        P1="trimmed/{sample}_1_paired.fastq",
        P2="trimmed/{sample}_2_paired.fastq",
        idx = config["genome_index_dir"]+"/genome.1.ht2"
    output:
        o1 = "final/aligned_{sample}.bam"
    params:
        index = config["genome_index_dir"] + "/genome"
    log:
        "log/rnaseq_{sample}.log"
    conda:
        "envs/rnaseq.yaml"
    shell:
        """
        set -x
        mkdir -p final
        hisat2 -p 8 --dta -x {params.index} -1 {input.P1} -2 {input.P2} | samtools sort -o {output.o1}
        rm -r {input.P1}
        rm -r {input.P2}
        """

rule use_featurecount:
    input:
        p1 = config["gtf_file"],
        p2 = "final/aligned_{sample}.bam"
    output:
        "final/{sample}_count.txt"
    log:
        "log/rnaseq_{sample}.log"
    conda:
        "envs/rnaseq.yaml"
    shell:
        """
        set -x
        featureCounts -S 2 -p -a {input.p1} -o {output} {input.p2}
        """