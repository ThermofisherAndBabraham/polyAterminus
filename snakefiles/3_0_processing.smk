
# ------------------------ Trimming ----------------------------------------- #

rule trim_adapters:
    input:
        tmp + "/{stem}_R1_001.fastq.gz",
        tmp + "/{stem}_R2_001.fastq.gz"
    output:
        LOGS + "/BBDUK/{stem}_contamination.log",
        temp(tmp + "/{stem}_R1_001Trimmed.fastq.gz"),
        temp(tmp + "/{stem}_R2_001Trimmed.fastq.gz")
    log:
        LOGS + "/BBDUK/trimming_{stem}.log"
    params:
        ref =       CONFIG["BBDUK"]["ref"],
        ktrim =     CONFIG["BBDUK"]["ktrim"],
        k =         CONFIG["BBDUK"]["k"],
        mink =      CONFIG["BBDUK"]["mink"],
        hdist =     CONFIG["BBDUK"]["hdist"],
        minlength = CONFIG["BBDUK"]["minlength"],
        qtrim =     CONFIG["BBDUK"]["qtrim"],
        trimq =     CONFIG["BBDUK"]["trimq"],
        add =       CONFIG["BBDUK"]["additional_params"],
        maxns =     CONFIG["BBDUK"]["maxns"],
        m =         MEMORY_JAVA
    threads:
        CONFIG["BBDUK"]["threads"]
    benchmark:
        BENCHMARKS + "/trimming_{stem}.log"
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


# ------------------------ Subsampling -------------------------------------- #

if CONFIG["SUBSAMPLING"]["subsample_to"]:

    rule subsample:
        input:
            tmp + "/{stem}_R1_001Trimmed.fastq.gz",
            tmp + "/{stem}_R2_001Trimmed.fastq.gz"
        output:
            temp(tmp + "/{stem}_R1_001subs.fastq.gz"),
            temp(tmp + "/{stem}_R2_001subs.fastq.gz")
        log:
            LOGS + "/SEQTK/R1/{stem}.log",
            LOGS + "/SEQTK/R2/{stem}.log"
        benchmark:
            BENCHMARKS + "/{stem}_subsampling.log"
        params:
            target_count = CONFIG["SUBSAMPLING"]["subsample_to"]
        threads:
            CONFIG["SUBSAMPLING"]["threads"]
        shell:
            "minimum={params.target_count}; " +
            "seqkit sample -s 11 -j {threads} --two-pass " +
            "-n {params.target_count} {input[0]} -o {output[0]} " +
            "2> {log[0]}; " +
            "seqkit sample -s 11 -j {threads} --two-pass " +
            "-n {params.target_count} {input[1]} -o {output[1]} " +
            "2> {log[1]}"

else:
    rule agregate_count:
        input:
            expand(tmp + "/{stem}_read_counts.txt", stem=STEMS)
        output:
            tmp + "/read_counts.txt"
        shell:
            "cat {input} > {output}"

    rule cp_read_counts:
        input:
            tmp + "/read_counts.txt"
        output:
            OUT + "/read_counts.txt"
        shell:
            "cp {input} {output}"

    rule count_reads:
        input:
            tmp + "/{stem}_R1_001.fastq.gz",
            tmp + "/{stem}_R1_001Trimmed.fastq.gz"
        output:
            tmp + "/{stem}_read_counts.txt"
        benchmark:
            BENCHMARKS + "/{stem}_count_reads.log"
        params:
            prefix = "{stem} ", nl = '''\n'''
        shell:
            "./scripts/fastq_num_reads.sh {input} > {output}; " +
            "sed -i -e 's/^/{params.prefix}/' {output} ; echo >>  {output}"

    rule subsample:
        input:
            tmp + "/{stem}_R1_001Trimmed.fastq.gz",
            tmp + "/{stem}_R2_001Trimmed.fastq.gz",
            OUT + "/read_counts.txt"
        output:
            temp(tmp + "/{stem}_R1_001subs.fastq.gz"),
            temp(tmp + "/{stem}_R2_001subs.fastq.gz")
        log:
            LOGS + "/SEQTK/R1/{stem}.log",
            LOGS + "/SEQTK/R2/{stem}.log"
        benchmark:
            BENCHMARKS + "/{stem}_subsampling.log"
        params:
            out = OUT
        threads:
            CONFIG["SUBSAMPLING"]["threads"]
        shell:
            """minimum=`python ./scripts/""" +
            """get_minimum_read_number.py --input {params.out}/read_counts.txt`;
            seqkit sample -s 11 -j {threads} --two-pass -n $minimum \
            {input[0]} -o {output[0]} 2> {log[0]};
            seqkit sample -s 11 -j {threads} --two-pass -n $minimum \
            {input[1]} -o {output[1]} 2> {log[1]};"""
