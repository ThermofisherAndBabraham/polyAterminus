# TrimmPolyA

This module provides methods for trimming polyA stretchers from reads wit the script mark_poly_A.jl


## Analysis

   1. Creation of encoded polyA database
      * Scan the reference genome (and transcriptome if annotation are available) searching for natural polyA stretches (min 20 bp)
      * Extract prefixes (by default min 10 bp) preceding the polyA and create FM-index for fast search against the prefix database.

   2. Detection of a polyA stretch
      * Scan from 3' end of R1 read and reverse complement of R2 read for a continues stretch of at least 10 bp
      of A.
      * A polyA stretch must occur to the second half  (towards 3' ) of the read (R1, R2 reverse complement) and the polyA stretches should be detected in both reads.

   3. The 3' of a polyA stretch having reads might have additional bases or even stretches of additional bases, examples:
      ```
      ...AAAAAAAAAAAAATAAAAAAAAAAAAAAAAACC
      ...AAAAAAAAAAAAATAAAATAAAAAAAAAAAAG
      ...AAAAAAAAAAAAATAAAATCCCCCCCCCGCCC
      ...AAAAAAAAAAAAAAGTTACTTTTTTTTTTTTTTT
      ```
   These bases are  “chopped off” the 3' of R1 read for futher processing
   4. Read trimming
      * Trim a polyA from the 3' of R1 counting A and allowing a certain number of mismatches per widow and overall per read (max 1 err per 10 bp window).
      * Check if a read doesn't' match a natural polyA stretches in a genome or transcripts.



## Result output

Results are stores in the four fastq.gz files:

   1. Reads without polyA stretches are outputted to two files with suffixes `..._R1_trimmedPolyA.fastq.gz` and `..._R2_trimmedPolyA.fastq.gz`A flag `-i` turn on inclusion of polyA trimmed reads.  

   2. Discarded polyA containing reads (ex. too short) reads (R1) are outputed to a files with suffix `..._discarded.fastq.gz`
   3. PolyA reads with trimmed polyA stretches are outputted to a file with suffix `_PolyA.fastq.gz`. The read names is modified indicating the trimmed length of a chopped off polyA stretch: `@<length of polyA>:A...`.  For example a read with trimmed off  22 bp poyA  stretch would have a name like this:
      * `@22:A:ST-E00243:412:HKCGMCCXY:6:2218:21958:59288`.  
