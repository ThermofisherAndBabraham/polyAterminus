
# ------------------------ Trimming ----------------------------------------- #

rule trim_adapters:
    input:
        tmp + "/{stem}_R1_001.fastq.gz",
        tmp + "/{stem}_R2_001.fastq.gz"
    output:
#        LOGS + "/BBDUK/{stem}_contamination.log",
        temp(tmp + "/{stem}_R1_001_val_1.fq.gz"),
        temp(tmp + "/{stem}_R2_001_val_2.fq.gz")
#    log:
#        LOGS + "/BBDUK/trimming_trim_galore_{stem}.log"
    params:
    threads:
        4
    benchmark:
        BENCHMARKS + "/trimming_trim_galore_{stem}.log"
    shell:
        "trim_galore --paired {input[0]} {input[1]} --output_dir {tmp}"

# removed the subsampling part

rule rename_trim_galore_files:
     input:
        tmp + "/{stem}_R1_001_val_1.fq.gz",
        tmp + "/{stem}_R2_001_val_2.fq.gz"
     output:
        tmp + "/{stem}_R1_001Trimmed.fastq.gz",
        tmp + "/{stem}_R2_001Trimmed.fastq.gz"
     shell:
        "cp -l {input[0]} {output[0]}; " +
        "cp -l {input[1]} {output[1]}"   

rule rename_to_subsample:
    input:
        tmp + "/{stem}_R1_001Trimmed.fastq.gz",
        tmp + "/{stem}_R2_001Trimmed.fastq.gz"
    output:
        tmp + "/{stem}_R1_001subs.fastq.gz",
        tmp + "/{stem}_R2_001subs.fastq.gz"
    benchmark:
        BENCHMARKS + "/{stem}_subsampling.log"
    shell:
        "cp -l {input[0]} {output[0]}; " +
        "cp -l {input[1]} {output[1]}" 