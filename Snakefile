#!python

import fnmatch
import yaml
import os.path
import sys
from shutil import copyfile
# --------------------------------- Config parsing ---------------------------- #

# Checks, if config file is specified through --configfile command. If not, uses
# config.yaml file.
if not config:
    with open("config.yaml") as stream:
        config = yaml.load(stream)

tmp_dir = config["tmp_dir"]
output_dir = config["output_dir"]

cmd = "if [ ! -d "+output_dir+" ]; then \
            mkdir "+output_dir+"; fi; \
       if [ ! -d "+tmp_dir+" ]; then \
            mkdir "+tmp_dir+"; fi;"
os.system(cmd)

if os.path.isdir(config["input_dir"]) == False:
    sys.exit("Provide input directory in config file.")

elif type(config["kmer_size"]) != int:
    sys.exit("Provide kmer size in config file.")

elif type(config["threads"]) != int:
    sys.exit("Provide threads number in config file.")

elif type(config["memory"]) != int:
    sys.exit("Provide memory size in config file.")

def get_stems(input_dir):
    stems = []
    for root, dirnames, filenames in os.walk(input_dir):
        for f in filenames:
            if f.find(".fastq.gz")>-1 or f.find(".fq.gz")>-1 : #and f.find("Z4-Leg-2_S8") > -1 :
                f = f.replace("_R1_001.fastq.gz", "")
                f = f.replace("_R2_001.fastq.gz", "")
                f = f.replace("1.fq.gz", "")
                f = f.replace("2.fq.gz", "")
                stems.append(f)
    return(list(set(stems)))

# Stems are taken from the source input directory

stems=get_stems(config["input_dir"])
print(stems)
max_threads = config["threads"]
max_picard_t = config["memory"]

# ------------------------------ Targets --------------------------- #

target = []


target += [",".join(expand("{tmp_dir}/{stem}_polyA_sorted.bam", \
                        stem=cls, tmp_dir=tmp_dir)) for cls in (stems)]
#                        stem=cls, tmp_dir=tmp_dir)) for cls in (stems)]
rule target:
    input: target

# ------------------------------ Get files -------------------------------- #
input_dir = config["input_dir"]

# ------------------------------ Anotation files -------------------------- #

config["annotation_files_prefix"]=annotation_files_prefix=os.path.splitext(config["reference_gtf"])[0]
config["genomic_dna_prefix"]=annotation_files_prefix=os.path.splitext(config["reference_dna"])[0]
config["genomic_dna_ercc_prefix"]=annotation_files_prefix=os.path.splitext(config["reference_dna"])[0]+"_ercc"
config["reference_dna_dir"]=annotation_files_prefix=os.path.dirname(config["reference_dna"])[0]
config["reference_dna_ercc"]=config["genomic_dna_prefix"]+"_ercc.fa"
config["reference_transcripts"]=config["genomic_dna_prefix"]+"_transcripts_ercc.fa"


#############################################################################
# ------------------------------ Analysis rules --------------------------- #
#############################################################################

# raw reads fastqc analysis

rule fastqc:
    input: input_dir+"/{stem}_R1_001.fastq.gz",
           input_dir+"/{stem}_R2_001.fastq.gz"
    output: "{tmp_dir}/{stem}_R1_001_fastqc.zip",
            "{tmp_dir}/{stem}_R2_001_fastqc.zip"
    params: kmer_size = config["kmer_size"]
    threads: 1
    wildcard_constraints:
        tmp_dir = config["tmp_dir"]
    shell: "fastqc {input} -o {tmp_dir} -t {threads} -k {params.kmer_size} " +
           ">> {output[0]} "

rule fastqc_post_proc:
    input: "{tmp_dir}/{stem}_R1_001Trimmed.fastq.gz",
           "{tmp_dir}/{stem}_R2_001Trimmed.fastq.gz"
    output: "{tmp_dir}/{stem}_R1_001Trimmed_fastqc.zip",
            "{tmp_dir}/{stem}_R2_001Trimmed_fastqc.zip"
    params: kmer_size = config["kmer_size"]
    threads: 1
    wildcard_constraints:
        tmp_dir = config["tmp_dir"]
    shell: "fastqc {input} -o {tmp_dir} -t {threads} -k {params.kmer_size} " +
           ">> {output[0]} "

rule multi_qc:
    input: expand("{tmp_dir}/{stem}_R1_001_fastqc.zip", stem=stems, tmp_dir=tmp_dir),
           expand("{tmp_dir}/{stem}_R2_001_fastqc.zip", stem=stems, tmp_dir=tmp_dir)
    output: "{output_dir}/fastqc_report_raw_reads.html",
            "{output_dir}/fastqc_report_raw_reads_data"
    shell: "multiqc {input} -o {output_dir} -n fastqc_report_raw_reads"

rule multi_qc_pre:
    input: expand("{tmp_dir}/{stem}_R1_001Trimmed_fastqc.zip", stem=stems, tmp_dir=tmp_dir),
           expand("{tmp_dir}/{stem}_R2_001Trimmed_fastqc.zip", stem=stems, tmp_dir=tmp_dir)
    output: "{output_dir}/fastqc_report_processed_reads.html",
            "{output_dir}/fastqc_report_processed_reads_data"
    shell: "multiqc {input} -o {output_dir} -n fastqc_report_processed_reads --force "


# ------------------------------ Trimming --------------------------- #

rule trim_adapters:
    input: input_dir + "/{stem}_R1_001.fastq.gz",
           input_dir + "/{stem}_R2_001.fastq.gz"
    output: temp("{tmp_dir}/{stem}_R1_001Trimmed.fastq.gz"),
            temp("{tmp_dir}/{stem}_R2_001Trimmed.fastq.gz"),
            "{tmp_dir}/{stem}_trim.settings"
    threads: config["threads"]
    params: bbmap_ref = config["bbmap_ref"]
    run:
            cmd = "nohup bbduk.sh in=%s in2=%s out=%s out2=%s " %(input[0],input[1],output[0],output[1])
            cmd += "ref=%s threads=%s " %(params.bbmap_ref,threads)
            cmd += "ktrim=r k=23 mink=11 hdist=1 minlength=100 maxns=1 tpe tbo   qtrim=r trimq=15 "
            cmd += ">%s; " %(output[2])    # Only forward reads are used vvv

            print(cmd)
            os.system(cmd)


# ------------------------------ Subsampling --------------------------- #

if config["subsample_to"] == None:
    rule subsample:
        input: "{tmp_dir}/{stem}_R1_001Trimmed.fastq.gz",
               "{tmp_dir}/{stem}_R2_001Trimmed.fastq.gz"
        output: "{tmp_dir}/{stem}_R1_001subs.fastq.gz",
                "{tmp_dir}/{stem}_R2_001subs.fastq.gz"
        params: tmp_dir="{tmp_dir}",
                f1="{stem}_R1_001subs.fastq.gz",
                f2="{stem}_R2_001subs.fastq.gz"
        threads: config["threads"]
        shell: "cd {params.tmp_dir} ; cp  ../{input[0]}  {params.f1} ; cp   ../{input[1]} {params.f2} ; cd .. "

else:
    rule subsample:
        input: "{tmp_dir}/{stem}_R1_001Trimmed.fastq.gz",
               "{tmp_dir}/{stem}_R2_001Trimmed.fastq.gz",
        output: "{tmp_dir}/{stem}_R1_001subs.fastq.gz",
                "{tmp_dir}/{stem}_R2_001subs.fastq.gz"
        params: target_count = config["subsample_to"]
        threads: config["threads"]
        shell: """minimum={params.target_count};
               seqkit sample -s 11 -j {threads} --two-pass \
               -n {params.target_count} {input[0]} -o {output[0]};
               seqkit sample -s 11 -j {threads} --two-pass \
               -n {params.target_count} {input[1]} -o {output[1]};"""


# ------------------------------ analysis of 3' coverages ----------------- #
#
rule  RSeQC_gene_body_coverage:
    input: bam = ["{tmp_dir}/{stem}_polyA_sorted.bam", \
                  "{tmp_dir}/{stem}_polyA_sorted.bam.bai"],
           bed = config["annotation_files_prefix"]+".bed"

    output: "{tmp_dir}/{stem}_gene_body_coverage.geneBodyCoverage.curves.png"
    params: prefix="{tmp_dir}/{stem}_gene_body_coverage"
    conda:
        "envs/rseqc.yaml"
    shell: "python2 ./scripts/geneBody_coverage_absolute.py -l 600 -i {input.bam[0]} -r {input.bed} -f png -o {params.prefix}"
#
# rule  RSeQC_read_distribution:
#     input: bam = ["{tmp_dir}/{stem}_polyA_sorted.bam", \
#                   "{tmp_dir}/{stem}_polyA_sorted.bam.bai"],
#            bed = config["annotation_files_prefix"]+".bed"
#
#     output: "{tmp_dir}/{stem}.read_distribution.txt"
#     params: sample = "{stem}"
#     shell: "read_distribution.py -i {input.bam[0]} -r {input.bed} " +
#            "| head -n 15  | tail -n 10 |  sed -e 's/^/{params.sample} /' " +
#            " >  {output[0]}"
#
rule  gtfToBed:
  input: config["annotation_files_prefix"]+".gtf"
  output: config["annotation_files_prefix"]+".bed"
  shell: " ./scripts/gtf2bed.pl {input} > {output}"
#
rule IndexPolyA:
    input:  "{tmp_dir}/{stem}_polyA_sorted.bam"
    output: "{tmp_dir}/{stem}_polyA_sorted.bam.bai"
    params:  threads = config["threads"]
    threads: max_threads
    shell: "sambamba index -p -t{params.threads} {input[0]}"
#
rule sortByNamePolyA:
    input: "{tmp_dir}/{stem}_polyA.bam"
    output: "{tmp_dir}/{stem}_polyA_sorted.bam"
    params: temp_dir_sambamba = config["scratch"],
            threads = config["threads"]
    threads: max_threads
    shell: "sambamba sort  --tmpdir={params.temp_dir_sambamba} -p " +
           "-t{params.threads} -o {output[0]} {input[0]}"

rule runSTARonPolyA:
    input: fastq="{tmp_dir}/{stem}_PolyA.fastq.gz",
               genomeDir=config["genomic_dna_ercc_prefix"]+"_STAR_annotation",
               genomeAnnotationFile=config["genomic_dna_ercc_prefix"] + \
                                "_STAR_annotation/Genome"
    threads: config["threads_star"]
    params: prefix = "{tmp_dir}/{stem}",
            starPrefix = "{tmp_dir}/{stem}Aligned.out",
            tmp = config["scratch"]+"/tmp_STAR_{stem}"
    output: "{tmp_dir}/{stem}_polyA.bam","{tmp_dir}/{stem}Log.final.out" #,"{tmp_dir}/{stem}_starTemp"
    run :
            cmd = "STAR --outTmpDir %s  --runThreadN %i " %(params.tmp, threads)
            cmd += "--genomeDir %s --runMode alignReads " %(input.genomeDir)
            cmd += "--outSAMunmapped Within --outFileNamePrefix %s " %(params.prefix)
            cmd += "--readFilesIn %s  --readFilesCommand zcat --genomeLoad NoSharedMemory " %(input.fastq)
            cmd += "--outSAMstrandField intronMotif  --outSAMattributes All  "
            cmd += " --outSAMtype BAM Unsorted --chimSegmentMin  20  " # --chimJunctionOverhangMin 20-outFilterMultimapNmax 20   --outFilterMismatchNmax 999  --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1 "
            cmd += " ; mv %s.bam %s " %(params.starPrefix, output[0])
            print(cmd)
            os.system(cmd)
#
rule trimPolyAReads:
        input:"{tmp_dir}/{stem}_R1_001subs.fastq.gz" , "{tmp_dir}/{stem}_R2_001subs.fastq.gz"
        output:"{tmp_dir}/{stem}_R1_trimmedPolyA.fastq.gz",
                "{tmp_dir}/{stem}_R2_trimmedPolyA.fastq.gz",
                "{tmp_dir}/{stem}_PolyA.fastq.gz",
                "{tmp_dir}/{stem}_discarded.fastq.gz"
        params: transcripts_fasta=config["transcripts_fasta"], output_stem="{tmp_dir}/{stem}"
        threads: config["threads_marker"]
        shell: "julia --depwarn=no scripts/mark_poly_A.jl  -i -p {threads} -a {input[0]} -b {input[1]} -o {params.output_stem} -c -r {params.transcripts_fasta} "

rule SearchA:
        input: "{tmp_dir}/{stem}_R1_001subs.fastq.gz",
        output: "{tmp_dir}/{stem}_R1_001subs_bbduk_selected_polyA.fastq",
                "{tmp_dir}/{stem}_polyAmasking.log"

        threads: config["threads"]
        params: bbmap_ref = config["polyA_ref"]
        run:
            cmd = " bbduk2.sh in=%s  outm=%s  " %(input[0],output[0])
            cmd += "fref=%s threads=%s rcomp=f " %(params.bbmap_ref,threads)
            cmd += " k=30  "
            cmd += " 2>&1 | tee  %s; " %(output[1])    # Only forward reads are used vvv
            print(cmd)
            os.system(cmd)


#
#
