
# ------------------------ Trimming polyA ----------------------------------- #
if (TRANSCRIPTS != None):
    rule PolyAAnalysis_trimm:
        input:
            PolyAAnalysis_trimm_input
        output:
            tmp + "/{stem}_R1_trimmedPolyA.fastq.gz",
            tmp + "/{stem}_R2_trimmedPolyA.fastq.gz",
            tmp + "/{stem}_PolyA.fastq.gz",
            tmp + "/{stem}_discarded.fastq.gz"
        benchmark:
            BENCHMARKS + "/{stem}_trim_polyA_reads.log"
        log:
            LOGS + "/TRIMMING-POLYA/{stem}.log"
        params:
            gz = TRANSCRIPTS,
            output_stem = tmp + "/{stem}"
        threads:
            JULIA_THREADS
        shell:
            "julia --depwarn=no PolyAAnalysis.jl/scripts/mark_poly_A.jl -i -p {threads} " +
            "-a {input[0]} -b {input[1]} -o {params.output_stem} " +
            "-r {params.gz} 2>&1 | tee -a {log}"

elif (gff != None and REFERENCE != None):
    rule PolyAAnalysis_trimm:
            input:
                PolyAAnalysis_trimm_input
            output:
                tmp + "/{stem}_R1_trimmedPolyA.fastq.gz",
                tmp + "/{stem}_R2_trimmedPolyA.fastq.gz",
                tmp + "/{stem}_PolyA.fastq.gz",
                tmp + "/{stem}_discarded.fastq.gz"
            benchmark:
                BENCHMARKS + "/{stem}_trim_polyA_reads.log"
            log:
                LOGS + "/TRIMMING-POLYA/{stem}.log"
            params:
                output_stem = tmp + "/{stem}",
                gff = gff,
                ref = REFERENCE
            threads:
                JULIA_THREADS
            shell:
                "julia --depwarn=no PolyAAnalysis.jl/scripts/mark_poly_A.jl -i -p {threads} " +
                "-a {input[0]} -b {input[1]} -o {params.output_stem} " +
                "-g {params.ref} -f {params.gff} 2>&1 | tee -a {log}"

elif (REFERENCE != None):
    rule PolyAAnalysis_trimm:
            input:
                PolyAAnalysis_trimm_input
            output:
                tmp + "/{stem}_R1_trimmedPolyA.fastq.gz",
                tmp + "/{stem}_R2_trimmedPolyA.fastq.gz",
                tmp + "/{stem}_PolyA.fastq.gz",
                tmp + "/{stem}_discarded.fastq.gz"
            benchmark:
                BENCHMARKS + "/{stem}_trim_polyA_reads.log"
            log:
                LOGS + "/TRIMMING-POLYA/{stem}.log"
            params:
                output_stem = tmp + "/{stem}",
                ref = REFERENCE
            threads:
                JULIA_THREADS
            shell:
                "julia --depwarn=no PolyAAnalysis.jl/scripts/mark_poly_A.jl -i -p {threads} " +
                "-a {input[0]} -b {input[1]} -o {params.output_stem} " +
                "-g {params.ref} 2>&1 | tee -a {log}"

else:
    wlogger.error(f'Reference , GFF3 or Transcripts file does not exist!')
    sys.exit("ERROR: reference, gff3 or transcripts files does not exist!")
