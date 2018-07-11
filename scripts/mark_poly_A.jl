#!/usr/bin/env julia

#= Custom library =#
using ArgParse.ArgParseSettings    #=                                    =#
using ArgParse.@add_arg_table      #= For parsing command-line options   =#
using ArgParse.parse_args          #=                                    =#
using CodecZlib
using IterTools
using JLD, HDF5
using FileIO

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
    use_cached_results - use already calculated results
"""
function get_polyA_prefixes_from_file(file::String;
    minimum_not_polyA::Int64=20,
    minimum_polyA_length::Int64=20,
    number_of_workers::Int64=4,
    use_cached_results::Bool=true)


    #file for caching
    jldFile=file*"_exracted_polyA_prefixes.jdl"

    if !isfile(jldFile) | !use_cached_results

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
        eval(macroexpand(quote @everywhere using PolyAAnalysis end))
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
        if use_cached_results

            save(File(format"JLD",jldFile), "all_result", all_result,compress=true)
            println(STDERR, "Colected unique polyA prefixes in known transcripts are saved in file $jldFile")
        end

    	close(file_stream)
    else
        println(STDERR ,"Loading colected unique polyA prefixes in known transcriptsfrom  file $jldFile")
        data = load(jldFile)
        all_result = data["all_result"]
        number_of_unque_prefixes=length(all_result)
        println(STDERR, "Loaded $number_of_unque_prefixes unique polyA prefixes in known transcripts")
    end

    return(all_result)
    println(STDERR,"colected prefixes are sto")
end

"""
finds and trims polyA having reads from a fastq pair
Arguments:
    fastq1 - file with forward reads in FASTQ format, can be gzipped
    fastq2 - file with reverse reads in FASTQ format, can be gzipped
    prefixes - array of prefixes of natural polyA
    number_of_workers - number of julia processes
    output_prefix - prefix for output files
    minimum_not_polyA - minimum length of not polyA strech in a read
    minimum_polyA_length - minimum length of polyA strech
    maximum_non_A_symbols - maximum numer of nonA symbols in polyA strech
    minimum_distance_from_non_poly - minimum distance of non A symbols in polyA strech from the nonA sequence
    maximum_distance_with_prefix_database - allowed Levenshtein distance between prefix of natural polyA in transcripts and the read
"""

function trim_polyA_from_files(
    fastq1::String,
    fastq2::String,
    prefixes::Array{String,1},
    number_of_workers::Int64,
    output_prefix::String,
    minimum_not_polyA::Int64,
    minimum_polyA_length::Int64,
    maximum_non_A_symbols::Int64,
    minimum_distance_from_non_poly_A::Int64,
    maximum_distance_with_prefix_database::Int64
    )

    #input streams

    file_stream1=open(fastq1,"r")
    if fastq1[length(fastq1)-2:end] == ".gz"
        file_stream1=GzipDecompressorStream(file_stream1)
    end
    file_stream2=open(fastq2,"r")
    if fastq2[length(fastq2)-2:end] == ".gz"
        file_stream2=GzipDecompressorStream(file_stream2)
    end
    #output streams
    fastqo1=GzipCompressorStream(open(output_prefix*"_R1_woPolyA.fastq.gz","w"))
    fastqo2=GzipCompressorStream(open(output_prefix*"_R2_woPolyA.fastq.gz","w"))
    fastqo_s=GzipCompressorStream(open(output_prefix*"_PolyA.fastq.gz","w"))
    fastqo_d=GzipCompressorStream(open(output_prefix*"_discarded.fastq.gz","w"))

    ostream1 = BioSequences.FASTQ.Writer(fastqo1,quality_header=false)
    ostream2 = BioSequences.FASTQ.Writer(fastqo2,quality_header=false)
    ostreams = BioSequences.FASTQ.Writer(fastqo_s,quality_header=false)
    ostreamd = BioSequences.FASTQ.Writer(fastqo_d,quality_header=false)


    # counter
    counter = convert(SharedArray, zeros(Int64, nworkers()))


    ct_pair=0
    ct_poly_a_proper=0
    ct_discarded=0


    for records in zip(collect(FASTQ.Reader(file_stream1)), collect(FASTQ.Reader(file_stream2)))

        revseq_rev_read=String(reverse_complement!(FASTQ.sequence(records[2])))
        seq_for_read=String(FASTQ.sequence(records[1]))
        ct_pair+=1
        #initial filtering
        polyA_detected=false
        if detect_polyA_in_a_string(seq_for_read,minimum_polyA_length,maximum_non_A_symbols) &&
            detect_polyA_in_a_string(revseq_rev_read,minimum_polyA_length,maximum_non_A_symbols)
                #tries sto esxtend polyA
                cosensus_extended_fwd_read=extend_poly_A(records[1],records[2])
                fqo_trimmed, has_proper_polyA, polyA_detected = trim_polyA_from_fastq_record(cosensus_extended_fwd_read,
                    minimum_not_polyA,
                    minimum_polyA_length,
                    maximum_non_A_symbols,
                    minimum_distance_from_non_poly_A)
                #check if the read is in pefixes list of naturall polyA sreches
                if has_proper_polyA
                    has_proper_polyA = check_polyA_prefixes(fqo_trimmed,prefixes,maximum_distance_with_prefix_database)
                end

        end
        #detect_polyA_in_fastq_record(records[1],minimum_polyA_length,maximum_non_A_symbols)
        if !polyA_detected
            write(fastqo1,records[1])
            write(fastqo2,records[2])
        else
            if has_proper_polyA
                write(fastqo_s,fqo_trimmed)
            else
                write(fastqo_d,cosensus_extended_fwd_read)
            end
        end

    end

    close(fastqo1)
    close(fastqo2)
    close(fastqo_s)
    close(fastqo_d)

end


function main(args)

    #= Command-line option parser =#
    arg_parse_settings = ArgParseSettings(description="Program trims polyA from the 3' end and modifies read name @[numberofAat3']_[originalname]")
    @add_arg_table arg_parse_settings begin
        "--output","-o"
            help="Output prefix"
            required = true
            arg_type = String

        "--fastq-f","-a"
            help="Input fatstq forward (R1) reads"
            required = true
            arg_type = String
        "--fastq-r","-b"
            help="Input fatstq reverse (R2) reads"
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

        "--use-precalculated-reference-transcripts-prefixes","-c"
            help="Load polyA prefixes from previous run"
            action = :store_true

        "--processes","-p"
            help="Number of additional julia workers for parallel procesing"
            required = true
            arg_type = Int64

    # testing julia --depwarn=no scripts/mark_poly_A.jl  -p -a  tmp/UHRR_100_Collibri_A_L001_R1_001subs.fastq.gz -b tmp/UHRR_100_Collibri_A_L001_R2_001subs.fastq.gz -o tmp/UHRR_100_Collibri_A_L001_R1_001subs_polyAmarked.fastq -r /media/root/ocean/Databases/hg38/gencode.v28.transcripts.fa.gz -c
     end

    parsed_args = parse_args(arg_parse_settings) #= In order to use variables, =#
                                                 #= use parsed_args["foo_bar"] =#
    #Parse reference file and collect prefixes of naturally occuring polyA prefixes


    polyA_prefixes=get_polyA_prefixes_from_file(parsed_args["reference-transcripts"],
    minimum_not_polyA=parsed_args["minimum-length"],
    minimum_polyA_length=parsed_args["minimum-polyA-length"],
    number_of_workers=parsed_args["processes"],
    use_cached_results=parsed_args["use-precalculated-reference-transcripts-prefixes"])
    trim_polyA_from_files(
        parsed_args["fastq-f"],
        parsed_args["fastq-r"],
        polyA_prefixes,
        parsed_args["processes"],
        parsed_args["output"],
        parsed_args["minimum-length"],
        parsed_args["minimum-polyA-length"],
        3, # maximum_non_A_symbols::Int64,
        1, # minimum_distance_from_non_poly::Int64
        3  # maximum_distance_with_prefix_database
        )


    #= Main code =#


end

main(ARGS)
