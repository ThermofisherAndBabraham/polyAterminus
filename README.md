## Description

Purpose of the work flow is to run RNA-seq trimming and alignment to a given reference steps and perform initial polyadenylation sites' analysis.

## Dependencies

* `julia v0.6.1`
* `snakemake 5.2.2`
* `BBMAP 38.22`
* `STAR 2.6.1a`
* `seqkit 0.8.1`
* `sambamba 0.6.6`


You can install  the dependencies manually or through conda environment as indicated below. If you choose to install the required software  manually please directly to the step 6 of the Workflow setup.

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
```bashr
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda config --add channels anaconda
    conda config --add channels r
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
   julia scripts/install_pkgs.jl
```




## Notes for `config.yaml`

* `Threads`: Consider is as how many parallel programs should be run on your PC.
* `memory`: restricts parallelism. Consider as threads.
* `memory_java`: for human 120 GB is recommended.
* `threads_julia`: usually no more than 6.
* `threads_star`:

## TODO

- [ ]  `output`: should be correct output for all necessary files.
- [ ]  `environment`: all applications should be added.
- [ ]  Tests for `MapPolyA.jl`
- [ ]  Check strands after parsing Bam in `MapPolyA.jl`
- [ ]  Test `rmdups()` in `MapPolyA.jl`
