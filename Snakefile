#!python

import fnmatch
import logging
import logging.config
import os.path
import re
import sys
import yaml
from shutil import copyfile


# ------------------------------ Include Snakefiles ------------------------- #

include: "./snakefiles/0_0_utilities.smk"
include: "./snakefiles/0_1_configuration.smk"
include: "./snakefiles/1_0_download.smk"
include: "./snakefiles/2_0_fastqc.smk"
include: "./snakefiles/3_0_processing.smk"
include: "./snakefiles/3_1_trimming_polyA.smk"
include: "./snakefiles/4_0_alignment.smk"
include: "./snakefiles/5_0_annotate_ts.smk"

init_log()
wlogger = logging.getLogger("custom_workflow")
wlogger.info(f'Input directory: [{CONFIG["INPUT-DIR"]}]')
wlogger.info(f'Using SAMPLES: {STEMS}')
wlogger.info(f'Lane regex for splitting filenames: [{CONFIG["LANE-REGEX"]}]')
wlogger.info(f'Tempory directory: [{CONFIG["TMP-DIR"]}]')
wlogger.info(f'Output direcoty: [{CONFIG["OUT-DIR"]}]')
wlogger.info(f'Reference: [{CONFIG["REFERENCE"]}]')
wlogger.info(f'GTF: [{CONFIG["GTF"]}]')
wlogger.info(f'GFF3: [{CONFIG["GFF3"]}]')

# ------------------------------ Targets ------------------------------------ #

rule all:
    input:
        expand("{d}/{s}_subSort.bam", d=ALIGN_DIR, s=STEMS),
        expand("{d}/STAR/{s}Log.final.out", d=OUT, s=STEMS),
        expand("{d}/ANNOTATE-POLYA/{s}_detected_polyA.tsv", d=OUT, s=STEMS),
        expand("{d}/ANNOTATE-POLYA/{s}_annotated_polyA.bed", d=OUT, s=STEMS),
        expand("{d}/ANNOTATE-POLYA/{s}_annotated_polyA_clusters_plus_coverage.bed", d=OUT, s=STEMS),
        MULTIQC_DIR + "/fastqc_report_raw_reads.html",
        MULTIQC_DIR + "/fastqc_report_trimmed_reads.html",
        MULTIQC_DIR + "/fastqc_report_processed_polyA_reads.html",


# ------------------------ Save run conf  ----------------------------------- #
rule git_log:
    output:
        LOGS + "/GIT-LOG/commit_used.log"
    shell:
        "echo $(git log) |  cut -d ' ' -f 2 &>> {output}"

rule write_config:
    output:
        LOGS + "/config.yaml"
    run:
        with open(logs + "/config.yaml", 'w') as outfile:
            yaml.dump(CONFIG, outfile, default_flow_style=False)
