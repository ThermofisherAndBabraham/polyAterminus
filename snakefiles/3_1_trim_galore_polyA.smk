
# ------------------------ Trimming polyA ----------------------------------- #

rule PolyAAnalysis_trimm:
    input:
        tmp + "/{stem}_R1_001Trimmed.fastq.gz",
        tmp + "/{stem}_R2_001Trimmed.fastq.gz"
    output:
        "{stem}_R1_001Trimmed_val_1.fq.gz",
        "{stem}_R2_001Trimmed_val_2.fq.gz"
    threads:
        JULIA_THREADS
    shell:
        "trim_galore --paired --polyA "+
        "{input[0]} {input[1]} "

rule rename_trimmed_files:
    input:
        "{stem}_R1_001Trimmed_val_1.fq.gz",
        "{stem}_R2_001Trimmed_val_2.fq.gz"
    output:
        tmp + "/{stem}_R1_trimmedPolyA.fastq.gz",
        tmp + "/{stem}_R2_trimmedPolyA.fastq.gz"
    shell:
        "mv {input[0]} {output[0]}; " +
        "mv {input[1]} {output[1]}"
        
rule filter_polyA:
    input:
        tmp + "/{stem}_R1_trimmedPolyA.fastq.gz"
    output:
        tmp + "/{stem}_PolyA.fastq.gz"
    shell:
        "zcat {input[0]} | grep -A {3} PolyA | grep -v ^-- > {output[0]}"
        