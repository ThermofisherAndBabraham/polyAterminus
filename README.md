## Description

Purpose of the work flow is to run complete primary and secondary analysis of RNA sequencing with Illumina machines and determine 3' ends of transcripts, PolyA length.

## Dependencies

* `julia v0.6`
* `snakemake 4.5.0`
* `BBMAP 36.92`
* `STAR 2.5.3a`
* `fastqc  v0.11.5`
* `multiqc 1.5`
* `AdapterRemoval 2.2.2`
* `seqkit 0.7.0`
* `sambamba 0.6.6`

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

## Install Julia Pkgs

* Run script: `julia scripts/install_pkgs.jl`
