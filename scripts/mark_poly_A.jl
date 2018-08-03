#!/usr/bin/env julia

#= Custom library =#
using ArgParse.ArgParseSettings    #=                                    =#
using ArgParse.@add_arg_table      #= For parsing command-line options   =#
using ArgParse.parse_args          #=                                    =#
using CodecZlib
using IterTools
using JLD, HDF5
using FileIO
using BioSequences
using BufferedStreams


push!(LOAD_PATH, ".")
using PolyAAnalysis





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
finds and trims polyA having reads from a fastq files pair
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
    minimum_poly_A_between - minimum length of polyA strech that might occure between non A symbols forming a fragment to be trimmed off (starting from the very 3' end)
    include_polyA - Includes output polyA sequences as pseudo pair end's  into output fastq

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
    maximum_distance_with_prefix_database::Int64,
    minimum_poly_A_between::Int64,
    include_polyA::Bool
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
    #output files

    fastqo1_f=output_prefix*"_R1_trimmedPolyA.fastq.gz"
    fastqo2_f=output_prefix*"_R2_trimmedPolyA.fastq.gz"
    fastqo_s_f=output_prefix*"_PolyA.fastq.gz"
    fastqo_d_f=output_prefix*"_discarded.fastq.gz"
    # @everywere test=open("test1.fq"

    fastqo1=open(fastqo1_f,"w")
    fastqo2=open(fastqo2_f,"w")
    fastqo_s=open(fastqo_s_f,"w")
    fastqo_d=open(fastqo_d_f,"w")

    #initiate and close
    flush(fastqo1)
    flush(fastqo2)
    flush(fastqo_s)
    flush(fastqo_d)
    close(fastqo1)
    close(fastqo2)
    close(fastqo_s)
    close(fastqo_d)

    "Opens output files in worker nodes"
    function open_output_files(name::String)
        global suffix = "_"*string(myid()-1)
        global n1=split(name,"\n")[1]
        eval(Main,:(fastqo1=GzipCompressorStream(open(n1*suffix,"a"))))
        global n2=split(name,"\n")[2]
        eval(Main,:(fastqo2=GzipCompressorStream(open(n2*suffix,"a"))))
        global n3=split(name,"\n")[3]
        eval(Main,:(fastqo_s=GzipCompressorStream(open(n3*suffix,"a"))))
        global n4=split(name,"\n")[4]
        eval(Main,:(fastqo_d=GzipCompressorStream(open(n4*suffix,"a"))))
    end

    "Closes output files in worker nodes"
    function close_output_files(name::String)
        eval(Main,:(close(fastqo1)))
        eval(Main,:(close(fastqo2)))
        eval(Main,:(close(fastqo_s)))
        eval(Main,:(close(fastqo_d)))
    end

    "Set  output file names"
    function set_output_files(name::String)
        global suffix = "_"*string(myid()-1)
        global n1=split(name,"\n")[1]
        eval(Main,:(fastqo1_f=n1*suffix))
        eval(Main,:(fastqo1_f_all=n1))
        global n2=split(name,"\n")[2]
        eval(Main,:(fastqo2_f=n2*suffix))
        eval(Main,:(fastqo2_f_all=n2))
        global n3=split(name,"\n")[3]
        eval(Main,:(fastqo_s_f=n3*suffix))
        eval(Main,:(fastqo_s_f_all=n3))
        global n4=split(name,"\n")[4]
        eval(Main,:(fastqo_d_f=n4*suffix))
        eval(Main,:(fastqo_d_f_all=n4))
    end

    "Merges file names"
    function merge_output_files(name::String)
        global suffix = "_"*string(myid()-1)
        global n1=split(name,"\n")[1]
        eval(Main,:(fastqo1_f=n1*suffix))
        eval(Main,:(fastqo1_f_all=n1))
        global n2=split(name,"\n")[2]
        eval(Main,:(fastqo2_f=n2*suffix))
        eval(Main,:(fastqo2_f_all=n2))
        global n3=split(name,"\n")[3]
        eval(Main,:(fastqo_s_f=n3*suffix))
        eval(Main,:(fastqo_s_f_all=n3))
        global n4=split(name,"\n")[4]
        eval(Main,:(fastqo_d_f=n4*suffix))
        eval(Main,:(fastqo_d_f_all=n4))

        fastqo1_f=Main.fastqo1_f
        fastqo1_f_all=Main.fastqo1_f_all
        fastqo2_f=Main.fastqo2_f
        fastqo2_f_all=Main.fastqo2_f_all
        fastqo_s_f=Main.fastqo_s_f
        fastqo_s_f_all=Main.fastqo_s_f_all
        fastqo_d_f=Main.fastqo_d_f
        fastqo_d_f_all=Main.fastqo_d_f_all

        ln=`cat $fastqo1_f >> $fastqo1_f_all ; rm $fastqo1_f`

        println("cat $fastqo1_f >> $fastqo1_f_all ; rm $fastqo1_f")
        run(ln)
        ln=`cat $fastqo2_f >> $fastqo2_f_all ; rm $fastqo2_f  `
        run(ln)
        ln=`cat $fastqo_s_f >> $fastqo_s_f_all; rm $fastqo_s_f `
        run(ln)
        ln=`cat $fastqo_d_f >> $fastqo_d_f_all ; rm $fastqo_d_f`
        run(ln)
    end



    # Indicate output files for worker nodes and open files
    Files_output = convert(SharedArray,collect(join([fastqo1_f,fastqo2_f,fastqo_s_f,fastqo_d_f],"\n")))
    @sync for (i, p) in enumerate(workers())
        @spawnat p open_output_files(String(Files_output))
    end

    @sync for (i, p) in enumerate(workers())
        @spawnat p open_output_files(String(Files_output))
    end








    debug_id="ST-E00243:413:HKCCHCCXY:7:2103:16234:11435"
    debug=false

    #Counters for specific pairs
    ct_all = convert(SharedArray, zeros(Int64, nworkers()))
    ct_pair_without_polyA = convert(SharedArray, zeros(Int64, nworkers()))
    ct_pair_with_proper_polyA = convert(SharedArray, zeros(Int64, nworkers()))
    ct_pair_with_discarded_polyA = convert(SharedArray, zeros(Int64, nworkers()))

    tic()
    @sync @parallel for records in collect(zip(collect(FASTQ.Reader(file_stream1)), collect(FASTQ.Reader(file_stream2)))) #@sync @parallel
        trim_polyA_from_fastq_pair(
                                    records[1],
                                    records[2],
                                    prefixes,
                                    minimum_not_polyA,
                                    minimum_polyA_length,
                                    maximum_non_A_symbols,
                                    maximum_distance_with_prefix_database,
                                    minimum_poly_A_between,
                                    include_polyA,
                                    ct_all,
                                    ct_pair_without_polyA,
                                    ct_pair_with_proper_polyA,
                                    ct_pair_with_discarded_polyA,
                                    debug=debug
                                    )
    end
    toc()
    # close output files for worker nodes
    Files_output = convert(SharedArray,collect(join([fastqo1_f,fastqo2_f,fastqo_s_f,fastqo_d_f],"\n")))
    @sync for (i, p) in enumerate(workers())
        @spawnat p close_output_files(String(Files_output))
    end


    pecentage_proper=round(100*sum(ct_pair_with_proper_polyA)/sum(ct_all),2)
    pecentage_discarded=round(100*sum(ct_pair_with_discarded_polyA)/sum(ct_all),2)
    println(STDERR,"Parsed reads: ",sum(ct_all))
    println(STDERR," % of having proper polyA: $pecentage_proper")
    println(STDERR," % of having discarded polyA: $pecentage_discarded")
    println(STDERR,"Merging files produced by polyA")
    Files_output = convert(SharedArray,collect(join([fastqo1_f,fastqo2_f,fastqo_s_f,fastqo_d_f],"\n")))
    @sync for (i, p) in enumerate(workers())
        @spawnat p merge_output_files(String(Files_output))
    end



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
            default = 10 #10
        "--reference-transcripts","-r"
            help="Reference transcripts in fasta format"
            required = true
            arg_type = String

        "--use-precalculated-reference-transcripts-prefixes","-c"
            help="Load polyA prefixes from previous run"
            action = :store_true
        "--incude-polyA-in-output","-i"
            help="Includes output polyA sequences as pseudo pair end's  into output fastq pair"
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

    addprocs(parsed_args["processes"])
    # Load modeules in the workers
    eval(macroexpand(quote @everywhere using BioSequences end))
    eval(macroexpand(quote @everywhere using CodecZlib end))
    eval(macroexpand(quote @everywhere using BufferedStreams end))
    eval(macroexpand(quote @everywhere push!(LOAD_PATH, ".") end))
    eval(macroexpand(quote @everywhere using PolyAAnalysis end))



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
        1, # maximum_non_A_symbols::Int64,
        1, # minimum_distance_from_non_poly::Int64
        3,  # maximum_distance_with_prefix_database
        4, # minimum_poly_A_between
        parsed_args["incude-polyA-in-output"]
        )






    #= Main code =#


end

main(ARGS)
