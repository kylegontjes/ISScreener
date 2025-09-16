# Author: Kyle Gontjes
# Date: 09-16-2025

configfile: "config/config.yaml"

import pandas as pd
import os
import re

PREFIX = config["prefix"]

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

if not os.path.exists("results/" + PREFIX):
    try:
        os.makedirs("results/" + PREFIX)
    except OSError as e:
        print(f"Error creating directory: {e}")

# Rule all information
rule all:
    input:
        final_bam = expand("results/{prefix}/{sample}/align_reads/{sample}_final.bam",prefix=PREFIX,sample=SAMPLE),
        ISFinder_out = expand("results/{prefix}/{sample}/panISa/{sample}_ISFinder.txt",prefix=PREFIX,sample=SAMPLE)

# Step 1: Trim raw fastq files using trimmomatic
rule trim_raw_reads:
    input:
        r1 = lambda wildcards: expand(str(config["input_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["input_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz"))
    output:
        r1_paired = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R1_trim_paired.fastq.gz",
        r2_paired = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R2_trim_paired.fastq.gz",
        r1_unpaired = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R1_trim_unpaired.fastq.gz",
        r2_unpaired = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R2_trim_unpaired.fastq.gz"
    params:
        num_threads=config["num_threads"],
        adapter_file_path=config["adapter_file_path"],
        seeds=config["seeds"],
        palindrome_clip_threshold=config["palindrome_clip_threshold"],
        simple_clip_threshold=config["simple_clip_threshold"],
        min_adapter_length=config["min_adapter_length"],
        keep_both_reads=config["keep_both_reads"],
        window_size=config["window_size"],
        window_size_quality=config["window_size_quality"],
        min_length=config["min_length"],
        head_crop_length=config["head_crop_length"]
    log:
        trim_log = "logs/{prefix}/{sample}/trimmomatic/{sample}_trimmomatic.log"
    singularity:
        "docker://staphb/trimmomatic:0.39"
    shell: 
        "trimmomatic PE {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} -threads {params.num_threads} ILLUMINACLIP:{params.adapter_file_path}:{params.seeds}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold}:{params.min_adapter_length}:{params.keep_both_reads} SLIDINGWINDOW:{params.window_size}:{params.window_size_quality}  MINLEN:{params.min_length} HEADCROP:{params.head_crop_length} &> {log.trim_log}"

# Step 2: Align reads using bwa and manipulate using samtools  
rule align_reads:
    input:
        r1_paired = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2_paired = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/{wildcards.sample}_R2_trim_paired.fastq.gz")
    output: 
        # Bam sorted files
        sorted_bam = temp(f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_sorted_aln.bam")
    params:
        reference_genome=config["reference_genome"],
        num_cores=config["num_cores"],
        outdir_temp = "results/{prefix}/{sample}/align_reads/{sample}_sorted_aln_temp"
    log: 
        sam_to_bam_log= "logs/{prefix}/{sample}/align_reads/{sample}_sam_to_bam.log"
    conda:
        "envs/bwa_samtools.yaml"
    shell:
        "bash ./bash_scripts/bwa.sh {input.r1_paired} {input.r2_paired} {params.reference_genome} {params.num_cores} {output.sorted_bam} {params.outdir_temp} {log.sam_to_bam_log} "
        
# Step 3: Remove duplicates and sort bam file w/ removed duplicates
rule remove_pcr_duplicates:
    input:
        sorted_bam = f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_sorted_aln.bam"
    output:
        # Bam dupliactes removed
        bam_duplicates_removed = f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_aln_marked.bam",
        # Picard duplicates_removed_name
        picard_metrics = f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_picard_metrics.txt",
    params:
        samtools_temp = "results/{prefix}/{sample}/align_reads/{sample}_final_temp"
    log:
        picard_log = "logs/{prefix}/{sample}/align_reads/{sample}_picard.log", 
    singularity:
        "docker://broadinstitute/picard:latest"
    shell:
        "java -jar /usr/picard/picard.jar MarkDuplicates -REMOVE_DUPLICATES true -INPUT {input.sorted_bam} -OUTPUT {output.bam_duplicates_removed} -METRICS_FILE {output.picard_metrics} -CREATE_INDEX false -VALIDATION_STRINGENCY LENIENT &> {log.picard_log}"

# Step 4: 
rule bam_sort:
    input:
        bam_duplicates_removed = f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_aln_marked.bam"
    output:
        # duplicates removed sorted bam
        final_bam = f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_final.bam"
    params:
        samtools_temp = "results/{prefix}/{sample}/align_reads/{sample}_final_temp"
    log:
        post_picard_sort = "logs/{prefix}/{sample}/align_reads/{sample}_post_picard_sort.log"
    singularity:
        "docker://staphb/samtools:1.19"
    shell:
        """
        samtools sort {input.bam_duplicates_removed} -m 500M -@ 0 -o {output.final_bam} -T {params.samtools_temp} &> {log.post_picard_sort}
        samtools index {output.final_bam}
        """

# Step 5: Run panISa & ISFinder
rule panISa:
    input:
        # sorted bam file
        final_bam = f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_final.bam"
    output:
        # panisa file
        panISa_out = f"results/{{prefix}}/{{sample}}/panISa/{{sample}}_panISa.txt",
        # isfinder file
        ISFinder_out = f"results/{{prefix}}/{{sample}}/panISa/{{sample}}_ISFinder.txt"
    log:
        panISA_log = "logs/{prefix}/{sample}/panISa/{sample}_panISa.log",
        ISFinder_log = "logs/{prefix}/{sample}/panISa/{sample}_ISFinder.log"
    conda:
        "envs/panISa.yaml"
    shell:
        """
        panISa.py {input.final_bam} -o {output.panISa_out} &> {log.panISA_log}
        ISFinder_search.py {output.panISa_out} -o {output.ISFinder_out} &> {log.ISFinder_log}
        """ 