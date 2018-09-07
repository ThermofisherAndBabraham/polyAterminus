#!/usr/bin/env julia

using ArgParse.ArgParseSettings
using ArgParse.@add_arg_table
using ArgParse.parse_args

push!(LOAD_PATH, ".")
using PolyAAnalysis

function main(args)
    arg_parse_settings = ArgParseSettings(description="Calculates uncovered windows number and GC bias.")
    @add_arg_table arg_parse_settings begin
        "--bam", "-b"
            nargs = '+'
            arg_type = String
            help = "Bam files separated by spaces."
            required = true
            dest_name = "bams"
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
    gffcollection = ParseGFF3(parsed_args["gff3file"])
    println("Parsed gff: ")
    toc()
    
    for bam in parsed_args["bams"]
        tic()
        treads, passreads, pareads, pclus = BamRead(bam)
        println("Parsed bam $bam: ")
        toc()
        tic()
        statdframe = PolACalculus(pclus)
        println("Gen dataframe: ")
        toc()
        statdframe[:Smaple] = basename(bam)
        statdframe[:TotalReads] = treads
        statdframe[:PassedReads] = passreads
        statdframe[:PolyAReads] = pareads
        statdframe = sort!(statdframe, [:Chrmosome, :Position])
        WrFrame(parsed_args["outfile"]*"_detected_polyA.tsv", statdframe, '\t')
        tic()
        intervalcolection = GetIntervalSet(statdframe)
        println("Got intervals in: ")
        toc()
        tic()
        joinedcollection = ItsectCollection(intervalcolection, gffcollection)
        println("Intersected in: ")
        toc()
        WrFrame(parsed_args["outfile"]*"_mapped_polyA.bed", joinedcollection, '\t')
    end
end

main(ARGS)
