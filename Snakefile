#!python

import fnmatch
import yaml
import os.path
import sys
from shutil import copyfile


# --------------------------------- Config parsing --------------------------- #

# Checks, if config file is specified through --configfile command. If not, uses
# config.yaml file.
if not config:
    with open("config.yaml") as stream:
        config = yaml.load(stream)

def get_stems(input_dir):
    stems = []
    for root, dirnames, filenames in os.walk(input_dir):
        for f in filenames:
            if f.find(".fastq.gz")>-1 or f.find(".fq.gz")>-1 :
                f = f.replace("_R1_001.fastq.gz", "")
                f = f.replace("_R2_001.fastq.gz", "")
                f = f.replace("1.fq.gz", "")
                f = f.replace("2.fq.gz", "")
                stems.append(f)
    return(list(set(stems)))

# Stems are taken from the source input directory
stems = get_stems(config["INPUT-DIR"])


# ------------------------------ Directories---------------------------------- #

if config["INPUT-DIR"][-1] == "/":
    input_dir = config["INPUT-DIR"][0:-2]
else:
    input_dir = config["INPUT-DIR"]

if config["TMP-DIR"][-1] == "/":
    tmp = config["TMP-DIR"][0:-2]
else:
    tmp = config["TMP-DIR"]

if config["OUT-DIR"][-1] == "/":
    out = config["OUT-DIR"][0:-2]
else:
    out = config["OUT-DIR"]

if config["SCRATCH"]:
    if config["SCRATCH"][-1] == "/":
        scratch = config["SCRATCH"][0:-2]
    else:
        scratch = config["SCRATCH"]
else:
    scratch = tmp

logs = out+"/"+"LOGS"


# ------------------------------ MACHINE ------------------------------------- #

threads = config["MACHINE"]["threads"]
star_threads = config["MAPPING"]["STAR"]["threads"]
julia_threads = config["MACHINE"]["threads_julia"]
memory_java = config["MACHINE"]["memory_java"]


# ------------------------------ REFERENCE ----------------------------------- #

reference = config["REFERENCE"]
gtf = config["GTF"]
gff = config["GFF3"]
transcripts = config["TRANSCRIPTS"]

# ------------------------------ Targets ------------------------------------- #

rule all:
    input:
        expand(out + "/STAR/{stem}_polyA.bam", out=out, stem=stems),
        expand(out + "/STAR/{stem}Log.final.out", out=out, stem=stems),
        expand(out + "/ANNOTATE-POLYA/{stem}_detected_polyA.tsv", out=out, stem=stems),
        expand(out + "/ANNOTATE-POLYA/{stem}_annotated_polyA.bed", out=out, stem=stems),
        out + "/fastqc_report_raw_reads.html",
        out + "/fastqc_report_processed_reads.html"

#############################################################################
# ------------------------------ Analysis rules --------------------------- #
#############################################################################

# ------------------------ Raw reads FASTQC  --------------------------------- #

rule fastqc:
    input:
        input_dir + "/{stem}_R1_001.fastq.gz",
        input_dir + "/{stem}_R2_001.fastq.gz"
    output:
        tmp + "/{stem}_R1_001_fastqc.zip",
        tmp + "/{stem}_R2_001_fastqc.zip"
    benchmark:
        "benchmarks/{stem}_raw_fastqc.log"
    params: kmer_size = config["FASTQC"]["kmer_size"]
    threads: 1
    wildcard_constraints:
        tmp = tmp
    shell:
        "fastqc {input} -o {tmp} -t {threads} -k {params.kmer_size} " +
        ">> {output[0]} "

rule fastqc_post_proc:
    input:
        tmp + "/{stem}_R1_001Trimmed.fastq.gz",
        tmp + "/{stem}_R2_001Trimmed.fastq.gz"
    output:
        tmp + "/{stem}_R1_001Trimmed_fastqc.zip",
        tmp + "/{stem}_R2_001Trimmed_fastqc.zip"
    benchmark:
        "benchmarks/{stem}_processed_fastqc.log"
    params:
        kmer_size = config["FASTQC"]["kmer_size"]
    threads: 1
    wildcard_constraints:
        tmp = tmp
    shell:
        "fastqc {input} -o {tmp} -t {threads} -k {params.kmer_size} " +
        ">> {output[0]} "

rule multi_qc:
    input:
        expand(tmp + "/{stem}_R1_001_fastqc.zip", stem=stems),
        expand(tmp + "/{stem}_R2_001_fastqc.zip", stem=stems)
    output:
        out + "/fastqc_report_raw_reads.html"
    shell:
        "multiqc {input} -o {out} -n fastqc_report_raw_reads --force"

rule multi_qc_pre:
    input:
        expand(tmp + "/{stem}_R1_001Trimmed_fastqc.zip", stem=stems),
        expand(tmp + "/{stem}_R2_001Trimmed_fastqc.zip", stem=stems)
    output:
        out + "/fastqc_report_processed_reads.html"
    shell:
        "multiqc {input} -o {out} -n fastqc_report_processed_reads --force "


# ------------------------ Trimming ------------------------------------------ #

rule trim_adapters:
    input:
        input_dir + "/{stem}_R1_001.fastq.gz",
        input_dir + "/{stem}_R2_001.fastq.gz"
    output:
        logs + "/BBDUK/{stem}_contamination.log",
        temp(tmp + "/{stem}_R1_001Trimmed.fastq.gz"),
        temp(tmp + "/{stem}_R2_001Trimmed.fastq.gz")
    log:
        logs + "/BBDUK/{stem}_trimming.log"
    params:
        ref =       config["BBDUK"]["ref"],
        ktrim =     config["BBDUK"]["ktrim"],
        k =         config["BBDUK"]["k"],
        mink =      config["BBDUK"]["mink"],
        hdist =     config["BBDUK"]["hdist"],
        minlength = config["BBDUK"]["minlength"],
        qtrim =     config["BBDUK"]["qtrim"],
        trimq =     config["BBDUK"]["trimq"],
        add =       config["BBDUK"]["additional_params"],
        maxns =     config["BBDUK"]["maxns"],
        m =         memory_java
    threads:
        config["BBDUK"]["threads"]
    benchmark:
        "benchmarks/trimming_{stem}.log"
    shell:
        "bbduk.sh in={input[0]} out={output[1]} " +
                "in2={input[1]} out2={output[2]} " +
                "ref={params.ref} threads={threads} " +
                "ktrim={params.ktrim} k={params.k} " +
                "mink={params.mink} hdist={params.hdist} " +
                "minlength={params.minlength} maxns={params.maxns} " +
                "qtrim={params.qtrim} trimq={params.trimq} " +
                "{params.add} threads={threads} " +
                "stats={output[0]} overwrite=t " +
                "-Xmx{params.m}g 2> {log}"


# ------------------------ Subsampling --------------------------------------- #

if config["SUBSAMPLING"]["run"]:
    if config["SUBSAMPLING"]["subsample_to"]:

        rule subsample:
            input:
                tmp + "/{stem}_R1_001Trimmed.fastq.gz",
                tmp + "/{stem}_R2_001Trimmed.fastq.gz"
            output:
                tmp + "/{stem}_R1_001subs.fastq.gz",
                tmp + "/{stem}_R2_001subs.fastq.gz"
            benchmark:
                "benchmarks/{stem}_subsampling.log"
            params:
                target_count = config["SUBSAMPLING"]["subsample_to"]
            threads:
                config["threads"]
            shell:
                "minimum={params.target_count}; " +
                "seqkit sample -s 11 -j {threads} --two-pass " +
                "-n {params.target_count} {input[0]} -o {output[0]}; " +
                "seqkit sample -s 11 -j {threads} --two-pass " +
                "-n {params.target_count} {input[1]} -o {output[1]}"

    else:
        rule agregate_count:
            input:
                expand(tmp + "/{stem}_read_counts.txt", stem=stems)
            output:
                tmp + "/read_counts.txt"
            shell:
                "cat {input} > {output}"

        rule cp_read_counts:
            input:
                tmp + "/read_counts.txt"
            output:
                out + "/read_counts.txt"
            shell:
                "cp {input} {output}"

        rule count_reads:
            input:
                input_dir + "/{stem}_R1_001.fastq.gz",
                tmp + "/{stem}_R1_001Trimmed.fastq.gz"
            output:
                tmp + "/{stem}_read_counts.txt"
            benchmark:
                "benchmarks/{stem}_count_reads.log"
            params:
                prefix = "{stem} ", nl = '''\n'''
            shell:
                "./scripts/fastq_num_reads.sh {input} > {output}; " +
                "sed -i -e 's/^/{params.prefix}/' {output} ; echo >>  {output}"

        rule subsample:
            input:
                tmp + "/{stem}_R1_001Trimmed.fastq.gz",
                tmp + "/{stem}_R2_001Trimmed.fastq.gz",
                out + "/read_counts.txt"
            output:
                temp(tmp + "/{stem}_R1_001subs.fastq.gz"),
                temp(tmp + "/{stem}_R2_001subs.fastq.gz")
            benchmark:
                "benchmarks/{stem}_subsampling.log"
            params:
                out = out
            shell:
                """minimum=`python ./scripts/""" +
                """get_minimum_read_number.py --input {params.out}/read_counts.txt`;
                seqkit sample -s 11 -j {threads} --two-pass -n $minimum \
                {input[0]} -o {output[0]};
                seqkit sample -s 11 -j {threads} --two-pass -n $minimum \
                {input[1]} -o {output[1]};"""

else:
    rule rename_to_subsample:
        input:
            tmp + "/{stem}_R1_001Trimmed.fastq.gz",
            tmp + "/{stem}_R2_001Trimmed.fastq.gz"
        output:
            tmp + "/{stem}_R1_001subs.fastq.gz",
            tmp + "/{stem}_R2_001subs.fastq.gz"
        benchmark:
            "benchmarks/{stem}_subsampling.log"
        threads:
            1
        shell:
            "cp {input[0]}  {output[0]}; " +
            "cp {input[1]} {output[1]}"


# ------------------------ analysis of 3' coverages -------------------------- #

rule IndexPolyA:
    input:
        tmp + "/{stem}_polyA_sorted.bam"
    output:
        tmp + "/{stem}_polyA_sorted.bam.bai"
    benchmark:
        "benchmarks/{stem}_IndexPolyA.log"
    params:
        threads = threads
    threads:
        threads
    shell:
        "sambamba index -p -t{params.threads} {input[0]}"

rule sortByNamePolyA:
    input:
        tmp + "/{stem}_polyA.bam"
    output:
        tmp + "/{stem}_polyA_sorted.bam"
    benchmark:
        "benchmarks/{stem}_sorByNamePolyA.log"
    params:
        scratch = scratch,
        threads = threads
    threads:
        threads
    shell:
        "sambamba sort  --tmpdir={params.scratch} -p " +
        "-t{params.threads} -o {output[0]} {input[0]}"

rule star_index:
    input:
        ref = reference,
        gtf = gtf
    output:
        directory(reference + "_STAR_annotation")
    benchmark:
        "benchmarks/star_index.log"
    threads:
        star_threads
    benchmark:
        "benchmarks/star_index.log"
    params:
        genomeDir = reference + "_STAR_annotation",
    shell:
        "if [ ! -d {params.genomeDir} ]; then " +
        " mkdir {params.genomeDir}; fi && " +
        " STAR --runMode genomeGenerate " +
            "--runThreadN {threads}  " +
            "--genomeDir {params.genomeDir} " +
            "--genomeFastaFiles {input.ref} " +
            "--sjdbGTFfile {input.gtf}"

rule alignment_polyA_reads:
    input:
        tmp + "/{stem}_R1_trimmedPolyA.fastq.gz",
        tmp + "/{stem}_R2_trimmedPolyA.fastq.gz",
        genomeDir = directory(reference + "_STAR_annotation")
    output:
        out + "/STAR/{stem}_polyA.bam",
        out + "/STAR/{stem}Log.final.out"
    benchmark:
        "benchmarks/{stem}_run_star.log"
    log:
        logs + "/STAR/{stem}Log.out"
    threads:
        star_threads
    params:
        prefix = tmp + "/{stem}",
        starPrefix = tmp + "/{stem}Aligned.out",
        tmp = scratch+"/tmp_STAR_{stem}",
        chimSegmentMin = config["MAPPING"]["STAR"]["chimSegmentMin"],
        additional = config["MAPPING"]["STAR"]["additional_params"]
    shell:
        "STAR --outTmpDir {params.tmp} --runThreadN {threads} " +
            "--genomeDir {input.genomeDir} --runMode alignReads " +
            "--outSAMunmapped Within " +
            "--outFileNamePrefix {params.prefix} " +
            "--readFilesIn {input[0]} {input[1]} --readFilesCommand zcat " +
            "--genomeLoad NoSharedMemory " +
            "--outSAMstrandField intronMotif --outSAMattributes All " +
            "--outSAMtype BAM Unsorted " +
            "--chimSegmentMin {params.chimSegmentMin} " +
            "{params.additional}; " +
            "mv {params.starPrefix}.bam {output[0]}; " +
            "mv {params.prefix}Log.final.out {output[1]}; " +
            "mv {params.prefix}Log.out {log}"

if (transcripts != None):
    rule trim_polyA_reads:
        input:
            tmp + "/{stem}_R1_001subs.fastq.gz" ,
            tmp + "/{stem}_R2_001subs.fastq.gz"
        output:
            tmp + "/{stem}_R1_trimmedPolyA.fastq.gz",
            tmp + "/{stem}_R2_trimmedPolyA.fastq.gz",
            tmp + "/{stem}_PolyA.fastq.gz",
            tmp + "/{stem}_discarded.fastq.gz"
        benchmark:
            "benchmarks/{stem}_trim_polyA_reads.log"
        log:
            logs + "/TRIMMING-POLYA/{stem}.log"
        params:
            gz = transcripts,
            output_stem = tmp + "/{stem}"
        threads:
            julia_threads
        shell:
            "julia --depwarn=no PolyAAnalysis.jl/scripts/mark_poly_A.jl -i -p {threads} " +
            "-a {input[0]} -b {input[1]} -o {params.output_stem} " +
            "-r {params.gz} 2>&1 | tee -a {log}"

elif (gff != None and reference != None):
    rule trim_polyA_reads:
            input:
                tmp + "/{stem}_R1_001subs.fastq.gz" ,
                tmp + "/{stem}_R2_001subs.fastq.gz"
            output:
                tmp + "/{stem}_R1_trimmedPolyA.fastq.gz",
                tmp + "/{stem}_R2_trimmedPolyA.fastq.gz",
                tmp + "/{stem}_PolyA.fastq.gz",
                tmp + "/{stem}_discarded.fastq.gz"
            benchmark:
                "benchmarks/{stem}_trim_polyA_reads.log"
            log:
                logs + "/TRIMMING-POLYA/{stem}.log"
            params:
                output_stem = tmp + "/{stem}",
                gff = gff,
                ref = reference
            threads:
                julia_threads
            shell:
                "julia --depwarn=no PolyAAnalysis.jl/scripts/mark_poly_A.jl -i -p {threads} " +
                "-a {input[0]} -b {input[1]} -o {params.output_stem} " +
                "-g {params.ref} -f {params.gff} 2>&1 | tee -a {log}"

elif (reference != None):
    rule trim_polyA_reads:
            input:
                tmp + "/{stem}_R1_001subs.fastq.gz" ,
                tmp + "/{stem}_R2_001subs.fastq.gz"
            output:
                tmp + "/{stem}_R1_trimmedPolyA.fastq.gz",
                tmp + "/{stem}_R2_trimmedPolyA.fastq.gz",
                tmp + "/{stem}_PolyA.fastq.gz",
                tmp + "/{stem}_discarded.fastq.gz"
            benchmark:
                "benchmarks/{stem}_trim_polyA_reads.log"
            log:
                logs + "/TRIMMING-POLYA/{stem}.log"
            params:
                output_stem = tmp + "/{stem}",
                ref = reference
            threads:
                julia_threads
            shell:
                "julia --depwarn=no PolyAAnalysis.jl/scripts/mark_poly_A.jl -i -p {threads} " +
                "-a {input[0]} -b {input[1]} -o {params.output_stem} " +
                "-g {params.ref} 2>&1 | tee -a {log}"

else:
    sys.exit("ERROR: reference, gff3 or transcripts files does not exist!")

rule annotate_polyA:
    input:
        out + "/STAR/{stem}_polyA.bam"
    output:
        out + "/ANNOTATE-POLYA/{stem}_detected_polyA.tsv",
        out + "/ANNOTATE-POLYA/{stem}_annotated_polyA.bed"
    benchmark:
        "benchmarks/{stem}_mapping_polyA.log"
    log:
        logs + "/ANNOTATE-POLYA/{stem}_mapped_polyA.log"
    threads:
        julia_threads
    params:
        gff = gff,
        pref = out + "/ANNOTATE-POLYA/{stem}",
        k = config["ANNOTATE-TS"]["k"],
        add_params = config["ANNOTATE-TS"]["additional_params"]
    shell:
        "export JULIA_NUM_THREADS={threads}; julia --depwarn=no " +
        "PolyAAnalysis.jl/scripts/annotate_polyA.jl -b {input} -o {params.pref} " +
        "-g {params.gff} -k {params.k} {params.add_params} &> {log}"
