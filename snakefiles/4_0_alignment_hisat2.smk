
# ------------------------ alignment ---------------------------------------- #

rule hisat_alignment:
    input:
        tmp + "/{stem}_R1_trimmedPolyA.fastq.gz",
        tmp + "/{stem}_R2_trimmedPolyA.fastq.gz",
    output:
        #temp(tmp + "/STAR/{stem}_polyA.sam"),
        tmp + "/STAR/{stem}_polyA.sam",
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
        additional = CONFIG["MAPPING"]["STAR"]["additional_params"],
        genome = CONFIG['HISAT2_GENOME'],
        splice_sites = CONFIG['HISAT2_SPLICES']
    shell:
        "hisat2 -p {threads} " +
            "-x {params.genome} --phred33-quals  " +
            "--no-mixed --no-discordant " +
            "--known-splicesite-infile {params.splice_sites} " +
            "-1 {input[0]} -2 {input[1]} " +
            "-S {output[0]} --summary-file {output[1]}"
   
# samtools view -bS -F 4 -F 8 -F 256 - > hisat_mapping/${input[0]}_hisat.bam "

rule samtobam:
    input:
        tmp + "/STAR/{stem}_polyA.sam"
    output:
        tmp + "/STAR/{stem}_polyA.bam"
    params:
        output_format = "bam"
    shell:
        "sambamba view --sam-input " + 
        "--format {params[0]} " +
        "--output-filename {output[0]} " +
        "{input[0]}"

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
