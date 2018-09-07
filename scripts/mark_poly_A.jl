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
function get_polyA_prefixes_from_file(file::Any, genomeFa::Any, gff::Any;
    minimum_not_polyA::Int64=20,
    minimum_polyA_length::Int64=20,
    number_of_workers::Int64=4,
    use_cached_results::Bool=true)

    #file for caching
    jldFile="exracted_polyA_prefixes.jdl"

    if !isfile(jldFile) | !use_cached_results
        # Open files and prepare decompression stream
        if (file != nothing)
            file_stream=open(file,"r")
        	if file[length(file)-2:end] == ".gz"
            	file_stream=GzipDecompressorStream(file_stream)
        	end
            collectedFasta = collect(FASTA.Reader(file_stream))
            close(file_stream)
        elseif (gff != nothing && genomeFa != nothing)
            collectedFasta = get_transcripts_from_gff(genomeFa, gff)
        else
            println(STDERR,"ERROR! Transcripts files")
        end
        println("Colleceted",collectedFasta)
        # Start julia worker processors

    	#Crate counter for progress nonitoring
        counter = convert(SharedArray, zeros(Int64, nworkers()))
    	# Arry to collect results for output
        all_result=Array{String,1}()
        #Parse transcripts
        time= @elapsed result= @parallel  (vcat) for record in collectedFasta
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
    chunk_fastq_analysis - Fastq pairs which are analysed together

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
    include_polyA::Bool;
    analysis_partition_size::Int64=100
    )

    #input streams



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

    #counter for jobs
    ct_jobs=0
    jub_submision_finished=false

    debug_id="ST-E00243:412:HKCGMCCXY:1:1120:14357:35432"
    debug=true

    #Counters for specific pairs
    ct_all = convert(SharedArray, zeros(Int64, nworkers()))
    ct_pair_without_polyA = convert(SharedArray, zeros(Int64, nworkers()))
    ct_pair_with_proper_polyA = convert(SharedArray, zeros(Int64, nworkers()))
    ct_pair_with_discarded_polyA = convert(SharedArray, zeros(Int64, nworkers()))


    const fastqpairs = RemoteChannel(()->Channel{ Array{Tuple{BioSequences.FASTQ.Record,BioSequences.FASTQ.Record},1} }(100));
    const fastqo_1_2_s_d=RemoteChannel(()->Channel{Tuple{String,String,String,String} }(100));

    @everywhere function trim_polyA_from_fastq_pair_pararell(
            fastq_pairs::RemoteChannel{Channel{Array{Tuple{BioSequences.FASTQ.Record,BioSequences.FASTQ.Record},1}}},
            fastqo_1_2_s_d::RemoteChannel{Channel{NTuple{4,String}}},
            ct_all::SharedArray{Int64,1},
            ct_pair_without_polyA::SharedArray{Int64,1},
            ct_pair_with_proper_polyA::SharedArray{Int64,1},
            ct_pair_with_discarded_polyA::SharedArray{Int64,1},
            prefixes::Array{String,1},
            minimum_not_polyA::Int64,
            minimum_polyA_length::Int64,
            maximum_non_A_symbols::Int64,
            maximum_distance_with_prefix_database::Int64,
            minimum_poly_A_between::Int64,
            include_polyA::Bool;
            debug_id="None",
            debug=false
            )
        while true
            fastq1_b=IOBuffer()
            fastq2_b=IOBuffer()
            fastq_s_b=IOBuffer()
            fastq_d_b=IOBuffer()
            for   (fastq1, fastq2) in take!(fastq_pairs)
                ct_all[(myid()-1)] += 1
                has_proper_polyA=false
                polyA_detected=false

                if debug_id == FASTQ.identifier(fastq1)
                    debug=true
                else
                    debug=false
                end



                fqo_trimmed, has_proper_polyA, polyA_detected = trim_polyA_from_fastq_pair(
                                                                fastq1,
                                                                fastq2,
                                                                prefixes,
                                                                minimum_not_polyA,
                                                                minimum_polyA_length,
                                                                maximum_non_A_symbols,
                                                                maximum_distance_with_prefix_database,
                                                                minimum_poly_A_between,
                                                                debug_id=debug_id,
                                                                debug=debug
                                                                )

                if debug
                    println("trim_polyA_from_fastq_pair output for $debug_id")
                    println("fqo_trimmed:\n $fqo_trimmed \n, has_proper_polyA: $has_proper_polyA, polyA_detected: $polyA_detected")
                end

                if !polyA_detected
            		ct_pair_without_polyA[(myid()-1)] += 1
                    println(fastq1_b,fastq1)
                    println(fastq2_b,fastq2)
                else
                    a=1
                    if has_proper_polyA
                        println(fastq_s_b,fqo_trimmed)
            			ct_pair_with_proper_polyA[(myid()-1)] += 1
                        #get reverse complement of the polyA read
                        if include_polyA
                            name=FASTQ.identifier(fqo_trimmed)
                            description=FASTQ.description(fqo_trimmed)
                            rev_quality=reverse(FASTQ.quality(fqo_trimmed))
                            rev_seq=BioSequences.reverse_complement!(FASTQ.sequence(fqo_trimmed))
                            fqo_trimmed_rev=FASTQ.Record(name, description, rev_seq, rev_quality)
                            println(fastq1_b,fqo_trimmed)
                            println(fastq2_b,fqo_trimmed_rev)
                        end
                    else

                        println(fastq_d_b,fastq1)
            			ct_pair_with_discarded_polyA[(myid()-1)] += 1
                    end
                end
            end

            #write out
            fastq1_b_s=String(fastq1_b)
            fastq2_b_s=String(fastq2_b)
            fastq_s_b_s=String(fastq_s_b)
            fastq_d_b_s=String(fastq_d_b)

            if fastq1_b_s != ""   fastq1_b_z=String(transcode(GzipCompressor,fastq1_b_s)) else fastq1_b_z="" end
            if fastq2_b_s != ""   fastq2_b_z=String(transcode(GzipCompressor,fastq2_b_s)) else fastq2_b_z="" end
            if fastq_s_b_s != ""   fastq_s_b_z=String(transcode(GzipCompressor,fastq_s_b_s)) else fastq_s_b_z="" end
            if fastq_d_b_s != ""   fastq_d_b_z=String(transcode(GzipCompressor,fastq_d_b_s)) else fastq_d_b_z="" end
            put!(fastqo_1_2_s_d,(fastq1_b_z,fastq2_b_z,fastq_s_b_z,fastq_d_b_z))
            close(fastq1_b)
            close(fastq2_b)
            close(fastq_s_b)
            close(fastq_d_b)
            end


    end





    function put_fastq_for_processing(fastq1::String,fastq2::String,analysis_partition_size::Int64)
        file_stream1=open(fastq1,"r")
        if fastq1[length(fastq1)-2:end] == ".gz"
            file_stream1=GzipDecompressorStream(file_stream1)
        end
        file_stream2=open(fastq2,"r")
        if fastq2[length(fastq2)-2:end] == ".gz"
            file_stream2=GzipDecompressorStream(file_stream2)
        end
        for records_pair in partition(zip(FASTQ.Reader(file_stream1),FASTQ.Reader(file_stream2)),analysis_partition_size)
            ct_jobs+=1 # global jobs counter
            #println("putted $ct_jobs ")
            put!(fastqpairs,collect(records_pair))
        end
        jub_submision_finished=true
    end

    @schedule put_fastq_for_processing(fastq1,fastq2,analysis_partition_size)






    for p in workers() # start tasks on the workers to process requests in parallel

              @async remote_do(trim_polyA_from_fastq_pair_pararell, p,
                                              fastqpairs,
                                              fastqo_1_2_s_d,
                                              ct_all,
                                              ct_pair_without_polyA,
                                              ct_pair_with_proper_polyA,
                                              ct_pair_with_discarded_polyA,
                                              prefixes,
                                              minimum_not_polyA,
                                              minimum_polyA_length,
                                              maximum_non_A_symbols,
                                              maximum_distance_with_prefix_database,
                                              minimum_poly_A_between,
                                              include_polyA;
                                              debug=debug,
                                              debug_id=debug_id
                                              )

    end
    ctout=0



    while  !((ctout == ct_jobs) & jub_submision_finished) || ((ctout==0) && (ct_jobs==0))# print out results
        ctout+=1
        fastqo1_s,fastqo2_s,fastqo_s_s,fastqo_d_s  = take!(fastqo_1_2_s_d)
        write(fastqo1,fastqo1_s)
        write(fastqo2,fastqo2_s)
        write(fastqo_s,fastqo_s_s)
        write(fastqo_d,fastqo_d_s)
        if mod(ctout*analysis_partition_size,1000)==0
            print(STDERR, "Parsed reads: ",sum(ct_all),
            " without polyA: ",sum(ct_pair_without_polyA),
            " with proper polyA: ",sum(ct_pair_with_proper_polyA),
            " with discarded polyA: ",sum(ct_pair_with_discarded_polyA),"\r")
        end
    end
    print(STDERR, "Parsed reads: ",sum(ct_all),
    " without polyA: ",sum(ct_pair_without_polyA),
    " with proper polyA: ",sum(ct_pair_with_proper_polyA),
    " with discarded polyA: ",sum(ct_pair_with_discarded_polyA),"\n")


    close(fastqo1)
    close(fastqo2)
    close(fastqo_s)
    close(fastqo_d)





    # function fill_with_fastq_pairs(n)
    #    for i in 1:n
    #        put!(jobs, i)
    #    end
    # end
    #
    #     # println(typeof(records_pairs))
    #     # exit()
    #     # println("iter within")
    #     # for p in records_pairs
    #     #     println("-----")
    #     #     println(p)
    #     #     println(typeof(p))
    #     # end
    #     # exit()
    #     #records_pairs SubArray{Tuple{BioSequences.FASTQ.Record,BioSequences.FASTQ.Record},1,Array{Tuple{BioSequences.FASTQ.Record,BioSequences.FASTQ.Record},1},Tuple{UnitRange{Int64}},true}
    #
    #     trim_polyA_from_fastq_pairs(
    #                                 records_pairs,
    #                                 prefixes,
    #                                 minimum_not_polyA,
    #                                 minimum_polyA_length,
    #                                 maximum_non_A_symbols,
    #                                 maximum_distance_with_prefix_database,
    #                                 minimum_poly_A_between,
    #                                 include_polyA,
    #                                 ct_all,
    #                                 ct_pair_without_polyA,
    #                                 ct_pair_with_proper_polyA,
    #                                 ct_pair_with_discarded_polyA,
    #                                 debug=debug
    #                                 )
    # end
    #
    # # close output files for worker nodes
    # Files_output = convert(SharedArray,collect(join([fastqo1_f,fastqo2_f,fastqo_s_f,fastqo_d_f],"\n")))
    # @sync for (i, p) in enumerate(workers())
    #     @spawnat p close_output_files(String(Files_output))
    # end
    #
    #
    # pecentage_proper=round(100*sum(ct_pair_with_proper_polyA)/sum(ct_all),2)
    # pecentage_discarded=round(100*sum(ct_pair_with_discarded_polyA)/sum(ct_all),2)
    # println(STDERR,"Parsed reads: ",sum(ct_all))
    # println(STDERR," % of having proper polyA: $pecentage_proper")
    # println(STDERR," % of having discarded polyA: $pecentage_discarded")
    # println(STDERR,"Merging files produced by polyA")
    #
    # for f in [fastqo1_f,fastqo2_f,fastqo_s_f,fastqo_d_f]
    #     cmd_strings=Array{String,1}()
    #     for (i, p) in enumerate(workers())
    #         nr_file=p-1
    #         fp="$f"*"_"*"$nr_file"
    #         #println(p)
    #         push!(cmd_strings," cat $fp >> $f ")
    #     end
    #     cmd_string=join(cmd_strings, " ; ")*" ; wait "
    #     run(`bash -c $cmd_string`)
    # end



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
            required = false
            arg_type = String
        "--reference-genome","-g"
            help="Reference genome in fasta format"
            required = false
            arg_type = String
        "--reference-gff","-f"
            help="Reference gff3 file"
            required = false
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


    println(parsed_args["reference-transcripts"], parsed_args["reference-genome"], parsed_args["reference-gff"])
    polyA_prefixes=get_polyA_prefixes_from_file(parsed_args["reference-transcripts"],
    parsed_args["reference-genome"],
    parsed_args["reference-gff"],
    minimum_not_polyA=parsed_args["minimum-length"],
    minimum_polyA_length=parsed_args["minimum-polyA-length"],
    number_of_workers=parsed_args["processes"],
    use_cached_results=parsed_args["use-precalculated-reference-transcripts-prefixes"])
    tic()
    # trim_polyA_from_files(
    #     parsed_args["fastq-f"],
    #     parsed_args["fastq-r"],
    #     polyA_prefixes,
    #     parsed_args["processes"],
    #     parsed_args["output"],
    #
    #     parsed_args["minimum-length"],
    #     parsed_args["minimum-polyA-length"],
    #     1, # maximum_non_A_symbols::Int64,
    #     1, # minimum_distance_from_non_poly::Int64
    #     3,  # maximum_distance_with_prefix_database
    #     4, # minimum_poly_A_between
    #     parsed_args["incude-polyA-in-output"]
    #     )
    toc()






    #= Main code =#


end

main(ARGS)
