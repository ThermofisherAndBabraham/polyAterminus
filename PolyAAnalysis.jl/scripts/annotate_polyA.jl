#!/usr/bin/env julia

using ArgParse.ArgParseSettings
using ArgParse.@add_arg_table
using ArgParse.parse_args

push!(LOAD_PATH, ".")
using PolyAAnalysis

function main(args)
    arg_parse_settings = ArgParseSettings(description = "Finds mapping positions of polyA sites, annotations and writes tsv and bed files.
                                                        Presumptions: 1) SE reads is used for alignment; 2) Reads were trimmed with 3EndD
                                                        PolyAAnalysis trimmer which adds tags to the read names.")
    @add_arg_table arg_parse_settings begin
        "--bam", "-b"
            arg_type = String
            help = "Bam file."
            required = true
            dest_name = "bam"
        "--out", "-o"
            arg_type = String
            help = "Prefix of output file."
            required = false
            dest_name = "outfile"
            default = "output"
        "--gff", "-g"
            arg_type = String
            help = "GFF3 file"
            required = true
            dest_name = "gff3file"
    end

    parsed_args = parse_args(arg_parse_settings)
    tic()
    gffcollection = parseGFF3(parsed_args["gff3file"])
    println("Parsed gff: ")
    toc()

    bam = parsed_args["bam"]
    tic()
    treads, passreads, pareads, pclus = readbam(bam)
    println("Parsed bam $bam: ")
    toc()
    tic()
    statdframe = stats_poly_a(pclus)
    println("Gen dataframe: ")
    toc()
    statdframe[:Smaple] = basename(bam)
    statdframe[:TotalReads] = treads
    statdframe[:PassedReads] = passreads
    statdframe[:PolyAReads] = pareads
    statdframe = sort!(statdframe, [:Chrmosome, :Position])
    wrframe(parsed_args["outfile"]*"_detected_polyA.tsv", statdframe, '\t')
    tic()
    intervalcolection = getintervals(statdframe)
    println("Got intervals in: ")
    toc()
    tic()
    joinedcollection = annotate_polya_sites(intervalcolection, gffcollection)
    println("Intersected in: ")
    toc()
    wrframe(parsed_args["outfile"]*"_annotated_polyA.bed", joinedcollection, '\t')

end

main(ARGS)
