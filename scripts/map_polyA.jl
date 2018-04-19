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

    gffcollection = ParseGFF3(parsed_args["gff3file"])

    for bam in parsed_args["bams"]
        treads, passreads, pareads, pclus = BamRead(bam)
        statdframe = PolACalculus(pclus)
        statdframe[:Smaple] = basename(bam)
        statdframe[:TotalReads] = treads
        statdframe[:PassedReads] = passreads
        statdframe[:PolyAReads] = pareads
        WrFrame(parsed_args["outfile"]*"detected_polyA.tsv", statdframe, '\t')
        intervalcolection = GetIntervalSet(statdframe)
        joinedcollection = ItsectCollection(intervalcolection, gffcollection)
        WrFrame(parsed_args["outfile"]*"_mapped_polyA.bed", joinedcollection, '\t')
    end
end

main(ARGS)
