# PolyATerminus
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square)](https://polya-terminus.s3-eu-central-1.amazonaws.com/index.html)

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
   ./miniconda/Miniconda3-latest-Linux-x86_64.sh &&
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

* Run `snakemake`:
```bash
    snakemake --configfile config.yaml -j 24 -k -p
```


## Notes for analysis configuration `config.yaml`

* `threads`: consider is as how many parallel programs should be running on your PC.
* `memory_java`: for human 60 GB is recommended.
* `threads_julia`: usually no more than 8.
* `threads_star`: max as you wish.
* `SCRATCH`: location of tmp directory for faster computation (SSD).
* `REFERENCE`: genome in fasta format.
* `GFF3` and `GTF`: should be the same version.

* If you want to normalize all samples by lowest read number found:
`SUBSAMPLING:`
    `run: true`
* If you want to subsample to specific read number:
`SUBSAMPLING:`
    `run: true`
    `subsample_to: N`
