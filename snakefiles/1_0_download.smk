
# ------------------------ Download  ---------------------------------------- #


rule clone_slimfastq:
    output:
        "dependencies/slimfastq/slimfastq"
    run:
        install_slimfastq()

rule generate_script_for_merging_files:
    input:
        "dependencies/slimfastq/slimfastq"
    output:
        tmp + "/{stem}R1_input_list.sh",
        tmp + "/{stem}R2_input_list.sh"
    run:
        get_sample_files_named_pipe_script(wildcards.stem)

rule get_fastq1:
    input:
        tmp + "/{stem}R1_input_list.sh"
    output:
        tmp + "/{stem}_R1_001.fastq.gz"
    threads: 1
    shell:
        "./{input} > {output}"

rule get_fastq2:
    input:
        tmp + "/{stem}R2_input_list.sh"
    output:
        tmp + "/{stem}_R2_001.fastq.gz"
    threads: 1
    shell:
        "./{input} > {output}"
