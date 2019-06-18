# PolyATerminus
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square)](https://polya-terminus.s3-eu-central-1.amazonaws.com/index.html)
https://travis-ci.com/ThermofisherAndBabraham/polyAterminus.svg?branch=master

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
    threads: consider it as how many parallel programs should be running on your PC.
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
```
