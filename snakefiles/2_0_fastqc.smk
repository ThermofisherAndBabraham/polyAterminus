
# ------------------------ Raw reads FASTQC  -------------------------------- #

rule fastqc:
    input:
        tmp + "/{stem}.fastq.gz"
    output:
        OUT + "/fastqc_report_htmls_zips/{stem}_fastqc.zip",
        OUT + "/fastqc_report_htmls_zips/{stem}_fastqc.html"
    log:
        LOGS + "/FASTQC/{stem}.log"
    benchmark:
        BENCHMARKS + "/fastqc_{stem}.log"
    params:
        kmer_size = CONFIG["FASTQC"]["kmer_size"],
        out_dir = OUT + "/fastqc_report_htmls_zips/"
    threads: 1
    shell:
        "fastqc {input} -o {params.out_dir} -t {threads} " +
        "-k {params.kmer_size} 2> {log}"

rule multiqc_pre:
    input:
        expand(OUT + "/fastqc_report_htmls_zips/{stem}_R{R}_001_fastqc.zip", stem=STEMS, R=['1','2'])
    output:
        MULTIQC_DIR + "/fastqc_report_raw_reads.html",
        directory(MULTIQC_DIR + "/fastqc_report_raw_reads_data")
    benchmark:
        BENCHMARKS + "/multi_qc.log"
    shell:
        "multiqc {input} -o "+MULTIQC_DIR+" -n fastqc_report_raw_reads"

rule multiqc_pro:
    input:
        expand(OUT + "/fastqc_report_htmls_zips/{stem}_R{R}_001Trimmed_fastqc.zip", stem=STEMS, R=['1','2'])
    output:
        MULTIQC_DIR + "/fastqc_report_trimmed_reads.html",
        directory(MULTIQC_DIR + "/fastqc_report_trimmed_reads_data")
    benchmark:
        BENCHMARKS + "/multi_qc.log"
    shell:
        "multiqc {input} -o "+MULTIQC_DIR+" -n fastqc_report_trimmed_reads"

rule multiqc_pro_polyA:
    input:
        expand(OUT + "/fastqc_report_htmls_zips/{stem}_R{R}_trimmedPolyA_fastqc.zip", stem=STEMS, R=['1','2'])
    output:
        MULTIQC_DIR + "/fastqc_report_processed_polyA_reads.html",
        directory(MULTIQC_DIR + "/fastqc_report_processed_polyA_reads_data")
    benchmark:
        BENCHMARKS + "/multi_qc_proc_polyA.log"
    shell:
        "multiqc {input} -o "+MULTIQC_DIR+" -n fastqc_report_processed_polyA_reads --force "
