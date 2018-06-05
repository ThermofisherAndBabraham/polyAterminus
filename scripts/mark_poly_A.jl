#!/usr/bin/env julia

#= Custom library =#
using ArgParse.ArgParseSettings    #=                                    =#
using ArgParse.@add_arg_table      #= For parsing command-line options   =#
using ArgParse.parse_args          #=                                    =#
using CodecZlib
using IterTools

push!(LOAD_PATH, ".")
using PolyAAnalysis
using  BioSequences



"""
returns minimum kmer streches of not polyA ins a supplied transcript sequences from a file
Arguments:
    sequence - string with polyA streches (type - BioSequences reference seq)
    minimum_not_polyA - minimum length of not polyA strech
    minimum_polyA_length - minimum length of polyA strech
    number_of_workers - number of julia processes


"""
function get_polyA_prefixes_from_file(file::String;
    minimum_not_polyA::Int64=20,
    minimum_polyA_length::Int64=20,
    number_of_workers::Int64=4)

	# Open files and prepare decompression stream
    file_stream=open(file,"r")
	if file[length(file)-2:end] == ".gz"
    	file_stream=GzipDecompressorStream(file_stream)
	end
    # Start julia worker processors
    addprocs(number_of_workers)
	# Load modeules in the workers
    eval(macroexpand(quote @everywhere using BioSequences end))
    eval(macroexpand(quote @everywhere push!(LOAD_PATH, ".") end))
    eval(macroexpand(quote @everywhere using PolyAanalysis end))
	#Crate counter for progress nonitoring
    counter = convert(SharedArray, zeros(Int64, nworkers()))
	# Arry to collect results for output
    all_result=Array{String,1}()
    #Parse transcripts
    time= @elapsed result= @parallel  (vcat) for record in collect(FASTA.Reader(file_stream))
        get_polyA_prefixes(record,minimum_not_polyA,
                          minimum_polyA_length,counter)
     end
	#get unique prefixes
    all_result=unique(result)
    number_of_unque_prefixes=length(all_result)
    println(STDERR, "Colected $number_of_unque_prefixes unique polyA prefixes in known transcripts in $time s")
    return(all_result)
	close(file_stream)
end


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
        "--minimum-polyA-length","-l"
            help="Minimum length of  polyA sequence"
            required = false
            arg_type = Int64
            default = 20
        "--reference-transcripts","-r"
            help="Reference transcripts in fasta format"
            required = true
            arg_type = String

        "--processes","-p"
            help="Number of additional julia workers for parallel procesing"
            required = true
            arg_type = Int64
    #polymorphism_limit::Float64=0.1, minimum_coverage::Int64=5 , mutations_number_limit
     end

    parsed_args = parse_args(arg_parse_settings) #= In order to use variables, =#
                                                 #= use parsed_args["foo_bar"] =#
    #Parse reference file and collect prefixes of naturally occuring polyA prefixes
    get_polyA_prefixes_from_file(parsed_args["reference-transcripts"],
    minimum_not_polyA=parsed_args["minimum-length"],
    minimum_polyA_length=parsed_args["minimum-polyA-length"],
    number_of_workers=parsed_args["processes"])
    exit()
    # println(STDERR, "Extracting prefixes of polyA streches from transcripts in file:\n$reference_transcripts")
    # polyA_in_transcripts = get_polyA_prefixes_from_file(
    # reference_transcripts,
    # minimum_not_polyA=parsed_args["minimum-length"],
    # minimum_polyA_length=50,
    # maximum_non_A_symbols=0,
    # minimum_distance_from_non_poly_A=10)
    # println(STDERR, @sprintf("Reference parsed, collected %i prefixes",length(polyA_in_transcripts)))

    #= Main code =#

    trim_polyA_file_records(parsed_args["output"],
    parsed_args["input"], parsed_args["minimum-length"])
end

main(ARGS)
