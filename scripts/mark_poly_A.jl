#!/usr/bin/env julia

#= Custom library =#
using ArgParse.ArgParseSettings    #=                                    =#
using ArgParse.@add_arg_table      #= For parsing command-line options   =#
using ArgParse.parse_args          #=                                    =#
using BioSequences

push!(LOAD_PATH, ".")
using PolyAanalysis




function main(args)

    #= Command-line option parser =#
    arg_parse_settings = ArgParseSettings(description="Program trims polyA from the 3' end and modifies read name @[numberofAat3']_[originalname]")
    @add_arg_table arg_parse_settings begin
        "--output","-o"
            help="Output fatstq"
            required = true
            arg_type = String
        "--input","-i"
            help="Input fatstq"
            required = true
            arg_type = String
        "--minimum-length","-m"
            help="Minimum length of not polyA sequence"
            required = false
            arg_type = Int64
            default = 20
    #polymorphism_limit::Float64=0.1, minimum_coverage::Int64=5 , mutations_number_limit
     end

    parsed_args = parse_args(arg_parse_settings) #= In order to use variables, =#
                                                 #= use parsed_args["foo_bar"] =#
    #= Main code =#
    trim_polyA_file_records(parsed_args["output"],
    parsed_args["input"],
    parsed_args["minimum-length"])
end

main(ARGS)
