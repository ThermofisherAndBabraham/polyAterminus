# PolyATerminus
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square)](https://polya-terminus.s3-eu-central-1.amazonaws.com/index.html)
[![Build Status](https://travis-ci.com/ThermofisherAndBabraham/polyAterminus.svg?branch=master)](https://travis-ci.com/ThermofisherAndBabraham/polyAterminus)

## Description
Purpose of the work flow is to run RNA-seq trimming and alignment to a given
reference steps and perform initial polyadenylation sites' analysis.

## Dependencies

* `julia v0.6.1`
* `snakemake 5.2.2`
* `BBMAP 38.22`
* `STAR 2.6.1a`
* `seqkit 0.8.1`
* `sambamba 0.6.6`


You can install  the dependencies manually or through conda environment as
indicated below. If you choose to install the required software  manually
please directly to the step 6 of the Workflow setup.

## Workflow setup

1. Install `conda`:
```bash
   wget -P miniconda https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh &&
   chmod 755 ./miniconda/Miniconda3-latest-Linux-x86_64.sh &&
   ./miniconda/Miniconda3-latest-Linux-x86_64.sh
```

2. Add path of `miniconda` to `.bashrc` if not selected option to add automatically during installation:
```bash
   cd ~/ && pth=$(pwd) &&
   echo "PATH=$PATH:$pth/miniconda3/bin" >> ~/.bashrc
```

3. Add channels to conda:
```bash
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda config --add channels anaconda
    conda config --add channels defaults
```

2. Create your `conda` environment:
 ```bash
    conda env create -f envs/3EndD.yaml -n 3EndD
 ```

5. Activate created environment:
```bash
    source activate 3EndD
```

6. Install the required julia packages:
```bash
    bash envs/build_julia_pkgs.sh
```

## Analysis

* Run `snakemake` as follows:
```bash
    snakemake --configfile  config.yaml -j 24 -k -p
```

This would ran an analysis on a small dataset matching human brain RNA-seq data of 21th chromosome.

## PolyAAnalysis tests

* Run `julia`:
```bash
    source activate 3EndD;
    julia PolyAAnalysis.jl/tests/runtests.jl
```

## Notes for analysis configuration `config.yaml`

* Configuration of your working environment:
```yaml
    INPUT-DIR: dir with fastq files
    SAMPLES: sample names, part of a file name. See example bellow.   
    LANE-REGEX: python regex for finding same sample files split by lanes. eg. use "L\\d\\d\\d_" to find any SAMPLE1_R[1,2]_L00[1,2,3]_001.fastq.gz SAMPLE2_R[1,2]_L00[1,2,3]_001.fastq.gz  in INPUT-DIR. Sample names should be unique.
    memory_java: for human 60 GB is recommended.
    threads_julia: specify threads for `julia`. Usually no more then 8. Some scripts scales only reading GZ | BAM files.
    threads_star: specify threads for `STAR` aligner.
    threads_sambamba: specify threads for sambamba tools.
    SCRATCH: location of tmp directory for faster computation (SSD).
    REFERENCE: genome in fasta format.
    GFF3: annotation files, should be the same version.
    GTF: annotation files, should be the same version.
    TRANSCRIPTS: Optional, generated from fasta and annotations if not provided.
```

* Supported file names:
```bash
    SAMPLE1_L001_R1_001.fastq.gz
    SAMPLE1_L001_R2_001.fastq.gz
```
```bash
    SAMPLE1_L001_1.fq.gz
    SAMPLE1_L001_2.fq.gz
```
```bash
    SAMPLE1_L001_R1_001.sfq
    SAMPLE1_L001_R2_001.sfq
```
```bash
    SAMPLE1_L001_1.sfq
    SAMPLE1_L001_1.sfq
```
SAMPLES and LANE-REGEX options in your config would be:
```yaml
    SAMPLES: SAMPLE1
    LANE-REGEX: "L\\d\\d\\d_"
```

* If you want to normalize all samples by lowest read number found:
```yaml
    SUBSAMPLING:
        run: true
```

* If you want to subsample to specific read number:
```yaml
    SUBSAMPLING:
        run: true
        subsample_to: N
```

* Additional parameters which are not listed in config can be passed as a string for BBDUK,
STAR and PolyA annotator programs.
```yaml
    additional_params: "passed as a string for specified programs."
```
* Specific parameters for polyA | termination sites annotator:
```yaml
    ANNOTATE-TS:
        k: distance from the cluster center allowed.
        m: minimum distance between clusters allowed. Two adjacent clusters with distance <= m will be merged.
        additional_params: pass -c|--cluster to run clustering, pass -v|--verbose to print proceeding of clustering.
        mappingquality: Only reads with greater than this mapping value will pass.
        strandedness:   List your SAMPLES bellow providing for each strandedness.
            HBR_100_Collibri_chr21small_A: If R1 read is the same as a gene sequence: "+". If R2 read is the same as a gene sequence: "-"
```

## Annotated BED files for termination sites

Workflow generates annotated transcription termination sites in `ANNOTATE-POLYA` directory.

For `*annotated_polyA_clusters.bed` files:

| Column Nr | Field Name         | Explanation |
| --------- | ------------------ | ----------- |
| 1         | Chr                | Chromosome  |
| 2         | Start              | Start (0-based) |
| 3         | End                | End |
| 4         | GeneName           | Gene Name |
| 5         | ClusterSize        | Detected number of reads representing this termination site. |
| 6         | Strand             | Strand |
| 7         | Feature            | Feature |
| 8         | ClusterCenter      | Position of cluster center. 1-based. |
| 9         | Biotype            | Biotype |
| 10        | ClusterMedian      | Median of frequencies of reads representing each termination site in cluster. |
| 11        | ClusterMean        | Mean of frequencies of reads representing each termination site in cluster. |
| 12        | ClusterMin         | Min of frequencies of reads representing each termination site in cluster. |
| 13        | ClusterMax         | Max of frequencies of reads representing each termination site in cluster. |
| 14        | Cluster1stQuartile | 1st quartile of frequencies of reads representing each termination site in cluster. |
| 15        | Cluster3rdQuartile | 3rd quartile of frequencies of reads representing each termination site in cluster. |
| 16        | TSMedian           | Median length polyA tail. |
| 17        | TSMean             | Mean length polyA tail. |
| 18        | TSMin              | Min length polyA tail. |
| 19        | TSMax              | Max length polyA tail. |
| 20        | TS1stQuartile      | 1st quartile of lengths of polyA tail. |
| 21        | TS3rdQuartile      | 3rd quartile of lengths of polyA tail. |

In addition for `*annotated_polyA_clusters_plus_coverage.bed` files

| Column Nr | Field Name         | Explanation |
| --------- | ------------------ | ----------- |
| 22        | - | Sambamba added Read count at Chr:Start-End |
| 23        | - | Sambamba added mean coverage at Chr:Start-End |
| 24        | - | Sambamba added Sample Name if any. |


For `*annotated_polyA.bed` files:

| Column Nr | Field Name         | Explanation |
| --------- | ------------------ | ----------- |
| 1         | Chr                | Chromosome |
| 2         | Start              | Start (0-based) |
| 3         | End                | End |
| 4         | GeneName           | Gene Name |
| 5         | Count              | Detected number of reads representing this termination site. |
| 6         | Strand             | Strand |
| 7         | Feature            | Feature |
| 8         | Median             | Median length polyA tail. |
| 9         | Min                | Min length polyA tail. |
| 10        | Max                | Max length polyA tail. |
| 11        | Biotype            | Biotype |
