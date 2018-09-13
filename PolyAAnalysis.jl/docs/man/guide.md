# Package Guide

## Installation

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

## Usage

PolyAAnalysis is designed to find and annotate polyA sites.

### Trimming and taging reads with polyA sequence.
1. If transcripts fasta are available trim and tag polyA reads in your fastq files by running:
```bash
    julia --depwarn=no PolyAAnalysis.jl/scripts/mark_poly_A.jl
        -i -p {threads}
        -a {input[0]} -b {input[1]}
        -o {output_prefix}
        -c -r {transcripts.fasta.gz}
```

2. Transcripts are not available, you can generate them from gff annotation
    and reference fasta and trim-tag polyA reads in your fastq files by running:
```bash
    julia --depwarn=no PolyAAnalysis.jl/scripts/mark_poly_A.jl
          -i -p {threads}
          -a {input[0]} -b {input[1]}
          -o {out_prefix}
          -c -g {reference.fasta}
          -f {reference.gff3}
```

3. If you have only fasta reference read tagging and trimming will be performed
    accordingly to reference sequences. Rich polyA sequences will not be trimmed
    from a read if relatively same sequence is in the genome.
```bash
    julia --depwarn=no PolyAAnalysis.jl/scripts/mark_poly_A.jl
          -i -p {threads}
          -a {input[0]} -b {input[1]}
          -o {output_prefix}
          -c -g {reference.fasta}
```

4. After trimming reads should be aligned to the reference using `STAR` for example.

### Annotating polyA sites.
1. To annotate detected polyA sites run:
```bash
    julia --depwarn=no PolyAAnalysis.jl/scripts/annotate_polyA.jl
          -b {input}
          -o {output_prefix}
          -g {reference.gff3}
```
All detected polyA sites with min, max, med are returned in tsv format and annotated bed file.
