
# ------------------------ alignment ---------------------------------------- #

rule star_index:
    input:
        ref = REFERENCE,
        gtf = GTF
    output:
        directory(REFERENCE + "_STAR_annotation")
    log:
        LOGS + "/STAR/INDEX.log"
    benchmark:
        BENCHMARKS + "/star_index.log"
    threads:
        STAR_THREADS
    benchmark:
        BENCHMARKS + "/star_index.log"
    params:
        genomeDir = REFERENCE + "_STAR_annotation",
    shell:
        "if [ ! -d {params.genomeDir} ]; then " +
        " mkdir {params.genomeDir}; fi && " +
        " STAR --runMode genomeGenerate " +
            "--runThreadN {threads}  " +
            "--genomeDir {params.genomeDir} " +
            "--genomeFastaFiles {input.ref} " +
            "--sjdbGTFfile {input.gtf} 2> {log}"

rule star_alignment:
    input:
        tmp + "/{stem}_R1_trimmedPolyA.fastq.gz",
        tmp + "/{stem}_R2_trimmedPolyA.fastq.gz",
        genomeDir = directory(REFERENCE + "_STAR_annotation")
    output:
        temp(tmp + "/STAR/{stem}_polyA.bam"),
        OUT + "/STAR/{stem}Log.final.out"
    benchmark:
        BENCHMARKS + "/{stem}_run_star.log"
    log:
        LOGS + "/STAR/{stem}Log.out"
    threads:
        STAR_THREADS
    params:
        prefix = tmp + "/{stem}",
        starPrefix = tmp + "/{stem}Aligned.out",
        tmp = scratch+"/tmp_STAR_{stem}",
        chimSegmentMin = CONFIG["MAPPING"]["STAR"]["chimSegmentMin"],
        additional = CONFIG["MAPPING"]["STAR"]["additional_params"]
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

rule sambamba_sort:
    input:
        tmp + "/STAR/{stem}_polyA.bam",
    output:
        ALIGN_DIR + "/{stem}_subSort.bam"
    log:
        LOGS + "/SAMBAMBA/SORT/{stem}.log"
    benchmark:
        BENCHMARKS + "/sort_bam_{stem}.log"
    params:
        temp_dir_sambamba = scratch
    threads:
        SAMBAMBA_THREADS
    shell:
        "sambamba sort --tmpdir={params.temp_dir_sambamba} -p " +
        "-t{threads}  -o {output[0]} {input[0]} 2> {log}"

rule sambamba_index:
    input:
        ALIGN_DIR + "/{stem}_subSort.bam"
    output:
        ALIGN_DIR + "/{stem}_subSort.bam.bai"
    log:
        LOGS + "/SAMBAMBA/INDEX/{stem}.log"
    benchmark:
        BENCHMARKS + "/sortIndex_{stem}.log"
    threads:
        SAMBAMBA_THREADS
    shell:
        "sambamba index -p -t{threads} {input[0]} 2> {log}"
