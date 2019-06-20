
# ------------------------ alignment ---------------------------------------- #

rule PolyAAnalysis_annotate:
    input:
        ALIGN_DIR + "/{stem}_subSort.bam"
    output:
        OUT + "/ANNOTATE-POLYA/{stem}_detected_polyA.tsv",
        OUT + "/ANNOTATE-POLYA/{stem}_annotated_polyA.bed",
        OUT + "/ANNOTATE-POLYA/{stem}_detected_clusters_polyA.tsv",
        OUT + "/ANNOTATE-POLYA/{stem}_annotated_polyA_clusters.bed"
    benchmark:
        BENCHMARKS + "/{stem}_mapping_polyA.log"
    log:
        LOGS + "/ANNOTATE-POLYA/{stem}_mapped_polyA.log"
    threads:
        CONFIG["ANNOTATE-TS"]["threads"]
    params:
        gff = gff,
        pref = OUT + "/ANNOTATE-POLYA/{stem}",
        k = CONFIG["ANNOTATE-TS"]["k"],
        m = CONFIG["ANNOTATE-TS"]["m"],
        add_params = CONFIG["ANNOTATE-TS"]["additional_params"]
    shell:
        "export JULIA_NUM_THREADS={threads}; julia --depwarn=no " +
        "PolyAAnalysis.jl/scripts/annotate_polyA.jl -b {input} -o {params.pref} " +
        "-g {params.gff} -k {params.k} -m {params.m} {params.add_params} &> {log}"

rule sambamba_depth_region:
    input:
        ALIGN_DIR + "/{stem}_subSort.bam",
        OUT + "/ANNOTATE-POLYA/{stem}_annotated_polyA_clusters.bed",
        ALIGN_DIR + "/{stem}_subSort.bam.bai"
    output:
        OUT + "/ANNOTATE-POLYA/{stem}_annotated_polyA_clusters_plus_coverage.bed"
    log:
        LOGS + "/SAMBAMBA/depth_{stem}.log"
    threads:
        SAMBAMBA_THREADS
    shell:
        "sambamba depth region -t {threads} -L {input[1]} " +
        "{input[0]} | sed '/^#/ d' > {output}"
