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
using FMIndexes
NEW_PATH=join(split(Base.source_path(),"/")[1:end-3],"/")
push!(LOAD_PATH, NEW_PATH)
using PolyAAnalysis


"""
get_polyA_prefixes_from_file(...)

Returns minimum kmer streches of not polyA ins a supplied transcript sequences from a file

Arguments:
    `file::String`: string with path to transcripts fasta file
    `genomeFa::String`: reference genome fasta file
    `gff::String`: reference gff file with transcripts and exons
    `minimum_not_polyA::Int64`: minimum length of not polyA strech
    `minimum_polyA_length::Int64`: minimum length of polyA strech
    `number_of_workers::Int64`: number of julia processes
    `use_cached_results::Bool`: use already calculated results
"""
function get_polyA_prefixes_from_file(file::Any, genomeFa::Any, gff::Any;
    minimum_not_polyA::Int64=20,
    minimum_polyA_length::Int64=20,
    number_of_workers::Int64=4,
    use_cached_results::Bool=true)::FMIndexes.FMIndex{7,}

    is_reference = true
    # file for caching
    if file != nothing
        jldFile = file * "_exracted_polyA_prefixes.jdl"
    elseif gff != nothing && genomeFa != nothing
        jldFile = gff * "_exracted_polyA_prefixes.jdl"
    elseif genomeFa != nothing
        jldFile = genomeFa * "_exracted_polyA_prefixes.jdl"
    else
        is_reference = false
    end


    if is_reference
        if (!isfile(jldFile) | !use_cached_results)
            # Open files and prepare decompression stream
            if file != nothing
                println(STDERR,"Reading transcripts")
                file_stream = open(file,"r")

            	if file[length(file)-2:end] == ".gz"
                	file_stream = GzipDecompressorStream(file_stream)
            	end
                collectedFasta = collect(FASTA.Reader(file_stream))
                close(file_stream)

            elseif gff != nothing && genomeFa != nothing
                collectedFasta = get_transcripts_from_gff(genomeFa, gff )
            elseif genomeFa != nothing
                println(STDERR,"Generating transcripts from genome")
                file_stream = open(genomeFa,"r")
                collectedFasta = collect(FASTA.Reader(file_stream))
                close(file_stream)
            else
                println(STDERR,"ERROR! Missing transcripts, gff, reference files")
            end
            # Start julia worker processors
        	# Crate counter for progress nonitoring
            counter = convert(SharedArray, zeros(Int64, nworkers()))
        	# Arry to collect results for output
            all_result = Array{String,1}()
            # Parse transcripts
            re = Regex("([ATGC]{$minimum_not_polyA})A{$minimum_polyA_length,}")
            println(STDERR,"Collecting PolyA")
            time = @elapsed result = @parallel  (vcat) for record in collectedFasta
                get_polyA_prefixes(record,minimum_not_polyA,
                                  minimum_polyA_length,re,counter)
            end
        	# get unique prefixes
            all_result = unique(result)
            number_of_unque_prefixes = length(all_result)
            println(STDERR, "Colected $number_of_unque_prefixes unique polyA prefixes in known transcripts in $time s")
            #prepare a concatenated string of the prefixes separated by : symbol
            concatenated_prefixes = join(all_result,':')
            println(length(concatenated_prefixes))
            index = FMIndex(concatenated_prefixes, r=4)

            if use_cached_results
                save(File(format"JLD",jldFile), "index", index,compress=true)
                println(STDERR, "FMI index of colected unique polyA prefixes in known transcripts are saved in file $jldFile")
            end

        else
            println(STDERR ,"Loading FMI index of colected unique polyA prefixes in known transcriptsfrom  file $jldFile")
            data = load(jldFile)
            index = data["index"]
            println(STDERR, "Loaded ....")
        end
    else
        println(STDERR,"Analysis will be done without checking for encoded polyA ")
        index = FMIndex(" ") #mark the situation with index of length 1
    end
    return index
end

function write_from_pipe(
    fastqo_1_2_s_d::RemoteChannel{Channel{NTuple{4,String}}},
    ct_output_chunks::SharedArray{Int64,1},
    ct_outputed_chunks::SharedArray{Int64,1},
    ct_submited::SharedArray{Int64,1},
    ct_finished::SharedArray{Int64,1},
    output_prefix::String)

    ctout = 0
    #output files
    fastqo1_f = output_prefix * "_R1_trimmedPolyA.fastq.gz"
    fastqo2_f = output_prefix * "_R2_trimmedPolyA.fastq.gz"
    fastqo_s_f = output_prefix * "_PolyA.fastq.gz"
    fastqo_d_f = output_prefix * "_discarded.fastq.gz"

    fastqo1 = open(fastqo1_f,"w")
    fastqo2 = open(fastqo2_f,"w")
    fastqo_s = open(fastqo_s_f,"w")
    fastqo_d = open(fastqo_d_f,"w")

    number_of_workers = length(ct_finished)

    while (ct_submited[1] == 0) || (ct_submited[1] > (sum(ct_outputed_chunks)+1))
        test = (ct_submited[1] == 0) || (ct_submited[1] > (sum(ct_outputed_chunks)+1))
        #println(ct_submited[1]," " ,sum(ct_outputed_chunks)+1," ",test)
        fastqo1_s,fastqo2_s,fastqo_s_s,fastqo_d_s = take!(fastqo_1_2_s_d)
        ctout += 1
        write(fastqo1,fastqo1_s)
        write(fastqo2,fastqo2_s)
        write(fastqo_s,fastqo_s_s)
        write(fastqo_d,fastqo_d_s)
        ct_outputed_chunks[1] = ctout -1 # one will be added after final write out


    end

    #println("__________________________________FINISHED", " ",sum(ct_output_chunks) ," ", ct_outputed_chunks[1])

    close(fastqo1)
    close(fastqo2)
    close(fastqo_s)
    close(fastqo_d)
    ct_outputed_chunks[1] += 1 #marks final write out
end

function monitor(
        ct_finished::SharedArray{Int64,1},
        ct_output_chunks::SharedArray{Int64,1},
        ct_outputed_chunks::SharedArray{Int64,1},
        ct_submited::SharedArray{Int64,1},
        ct_pair_without_polyA::SharedArray{Int64,1},
        ct_pair_with_proper_polyA::SharedArray{Int64,1},
        ct_pair_with_discarded_polyA::SharedArray{Int64,1},
        final_monitoring_message::SharedArray{Int64,1},
        ct_all::SharedArray{Int64,1})
    number_of_workers = length(ct_finished)
    sleep_time = 5
    sleep_ct = 0
    finished = 0
    finished_first = 0
    end_time = now()
    start_time = now()

    while  (ct_submited[1]  < sum(ct_outputed_chunks)) || (sum(ct_outputed_chunks) == 0)
            sleep(sleep_time)
            old_finished = finished
            finished = sum(ct_all)
            end_time = now()
            sleep_ct += 1

            if sleep_ct == 1
                start_time = end_time
                finished_first = finished
            end

            speed = round(((finished-old_finished)/sleep_time),2)
            print(STDERR, "Currently parsed reads: ",finished," ; ","at speed $speed read pairs/s ; ",
            " without polyA: ",sum(ct_pair_without_polyA)," ; ",
            " with proper polyA: ",sum(ct_pair_with_proper_polyA)," ; ",
            " with discarded polyA: ",sum(ct_pair_with_discarded_polyA)," ; ",
            " output chunks: ",sum(ct_outputed_chunks),"\r")
    end


    elapsed_time = convert(Int64,Dates.value(end_time-start_time))
    speed = round(((finished-finished_first)/elapsed_time*1000),2)
    print(STDOUT, "In total parsed reads: ",finished," ; ","at average speed $speed read pairs/s ; ",
    " without polyA: ",sum(ct_pair_without_polyA)," ; ",
    " with proper polyA: ",sum(ct_pair_with_proper_polyA)," ; ",
    " with discarded polyA: ",sum(ct_pair_with_discarded_polyA),"                      \n")
    final_monitoring_message[1] = 1
end

"""
trim_polyA_from_files(...)

finds and trims polyA having reads from a fastq files pair

Arguments:
    `fastq1::String`: file with forward reads in FASTQ format, can be gzipped
    `fastq2::String`: file with reverse reads in FASTQ format, can be gzipped
    `prefixes::FMIndexes`: array of prefixes of natural polyA
    `number_of_workers::Int64`: number of julia processes
    `output_prefix`: prefix for output files
    `minimum_not_polyA::Int64`: minimum length of not polyA strech in a read
    `minimum_polyA_length::Int64`: minimum length of polyA strech
    `maximum_non_A_symbols::Int64`: maximum numer of nonA symbols in polyA strech
    `maximum_distance_with_prefix_database::Int64`: allowed Levenshtein distance between prefix of natural polyA in transcripts and the read
    `minimum_poly_A_between::Int64`: minimum length of polyA strech that might occure between non A symbols forming a fragment to be trimmed off (starting from the very 3' end)
    `include_polyA::String`: includes output polyA sequences as pseudo pair end's  into output fastq
"""
function trim_polyA_from_files(
    fastq1::String,
    fastq2::String,
    prefixes::FMIndexes.FMIndex{7,},
    number_of_workers::Int64,
    output_prefix::String,
    minimum_not_polyA::Int64,
    minimum_polyA_length::Int64,
    maximum_non_A_symbols::Int64,
    maximum_distance_with_prefix_database::Int64,
    minimum_poly_A_between::Int64,
    include_polyA::Bool;
    )

    # Setup exraction named pipe
    r1_filename = basename(fastq1)
    r2_filename = basename(fastq2)
    rand_suffix = randstring()
    r1_filename_tmp_prefix_ini = "/tmp/" * r1_filename * rand_suffix
    r2_filename_tmp_prefix_ini = "/tmp/" * r2_filename * rand_suffix
    cmd1 = "mkfifo $r1_filename_tmp_prefix_ini ; zcat  $fastq1 > $r1_filename_tmp_prefix_ini & "
    cmd2 = "mkfifo $r2_filename_tmp_prefix_ini ; zcat  $fastq2 > $r2_filename_tmp_prefix_ini & "
    run(`bash -c $cmd1`)
    run(`bash -c $cmd2`)
    f1_file = open(r1_filename_tmp_prefix_ini,"r")
    f2_file = open(r2_filename_tmp_prefix_ini,"r")

    #prepare counters for monitoring of execution
    #counter for jobs
    ct_jobs = 0
    debug_id = "ST-E00243:413:HKCCHCCXY:8:2224:31791:40600"
    debug = true

    #Counters for specific pairs and run monitors
    ct_all = convert(SharedArray, zeros(Int64, nworkers()))
    ct_pair_without_polyA = convert(SharedArray, zeros(Int64, nworkers()))
    ct_pair_with_proper_polyA = convert(SharedArray, zeros(Int64, nworkers()))
    ct_pair_with_discarded_polyA = convert(SharedArray, zeros(Int64, nworkers()))
    ct_output_chunks = convert(SharedArray, zeros(Int64, nworkers()))
    ct_outputed_chunks = convert(SharedArray, zeros(Int64, 1)) #collect number of outputed pairs
    # a special counter to check that a worker has fnished processing its chink of data.
    # value 1 marks that it has finished
    ct_finished = convert(SharedArray, zeros(Int64, nworkers()))
    # marks finishing of submision for analysis processes
    ct_submited =  convert(SharedArray, zeros(Int64, 1))
    #checks if final monitoring mesage is printed out
    final_monitoring_message = convert(SharedArray, zeros(Int64, 1))

    min_l = Int64(minimum_polyA_length / 2)
    re = Regex("(A+[GTC])?(A{$min_l,})([GTC]A+)?")

    const fastqo_1_2_s_d = RemoteChannel(()->Channel{Tuple{String,String,String,String} }(100));



    @schedule write_from_pipe(
        fastqo_1_2_s_d,
        ct_output_chunks,
        ct_outputed_chunks,
        ct_submited,
        ct_finished,
        output_prefix)

    @schedule monitor(
        ct_finished,
        ct_output_chunks,
        ct_outputed_chunks,
        ct_submited,
        ct_pair_without_polyA,
        ct_pair_with_proper_polyA,
        ct_pair_with_discarded_polyA,
        final_monitoring_message,
        ct_all
    )

    #setup run monitoring

    println(STDERR,"Prefix $output_prefix will be used for output")

    chunk_size = 10000
    ct_partition = 0
    ct_reads = 0
    chunk = Array{Tuple{NTuple{4,String},NTuple{4,String}},1}()
    ct_chunks = 0
    worker_nr=1

    for part in zip(partition(eachline(f1_file, chomp = false),4),partition(eachline(f2_file, chomp=false),4))
        push!(chunk,part)
        ct_partition += 1
        ct_reads += 1
        if ct_partition == chunk_size
             ct_chunks  += 1
             remote_do(trim_polyA_from_fastq_pair_pararell, worker_nr+1,
                                            chunk,
                                            fastqo_1_2_s_d,
                                            ct_all,
                                            ct_pair_without_polyA,
                                            ct_pair_with_proper_polyA,
                                            ct_pair_with_discarded_polyA,
                                            ct_finished,
                                            ct_output_chunks,
                                            prefixes,
                                            minimum_not_polyA,
                                            minimum_polyA_length,
                                            maximum_non_A_symbols,
                                            maximum_distance_with_prefix_database,
                                            minimum_poly_A_between,
                                            include_polyA,
                                            re;
                                            debug=debug,
                                            debug_id=debug_id
                                            )
            worker_nr += 1
            if worker_nr > number_of_workers
                worker_nr = 1
            end
            chunk = Array{Tuple{NTuple{4,String},NTuple{4,String}},1}()
            ct_partition = 0
        end

    end
    # final chunk analysis
    ct_chunks  += 1
    remote_do(trim_polyA_from_fastq_pair_pararell, worker_nr+1,
                                    chunk,
                                    fastqo_1_2_s_d,
                                    ct_all,
                                    ct_pair_without_polyA,
                                    ct_pair_with_proper_polyA,
                                    ct_pair_with_discarded_polyA,
                                    ct_finished,
                                    ct_output_chunks,
                                    prefixes,
                                    minimum_not_polyA,
                                    minimum_polyA_length,
                                    maximum_non_A_symbols,
                                    maximum_distance_with_prefix_database,
                                    minimum_poly_A_between,
                                    include_polyA,
                                    re;
                                    debug=debug,
                                    debug_id=debug_id
                                    )
    println(STDERR, "Submited for analysis $ct_reads reads in $ct_chunks chunks")
    ct_submited[1] = ct_chunks
    #wait till all reads will be outputed or final mesage of monitoring task will be outputed
    while (ct_submited[1] == 0) || (ct_submited[1] > (sum(ct_outputed_chunks))) || (final_monitoring_message[1]==0)
        sleep(1)
    end
    println("Analysis finished")

end

function main(args)

    #= Command-line option parser =#
    arg_parse_settings = ArgParseSettings(description="Program trims polyA from the 3' end and modifies read name @[numberofAat3']_[originalname]")
    @add_arg_table arg_parse_settings begin
        "--output","-o"
            help = "Output prefix"
            required = true
            arg_type = String
        "--fastq-f","-a"
            help = "Input fatstq forward (R1) reads"
            required = true
            arg_type = String
        "--fastq-r","-b"
            help = "Input fatstq reverse (R2) reads"
            required = true
            arg_type = String
        "--minimum-length","-m"
            help = "Minimum length of not polyA sequence"
            required = false
            arg_type = Int64
            default = 20
        "--minimum-polyA-length","-l"
            help = "Minimum length of  polyA sequence"
            required = false
            arg_type = Int64
            default = 10 #10
        "--reference-transcripts","-r"
            help = "Reference transcripts in fasta format"
            required = false
            arg_type = String
        "--reference-genome","-g"
            help = "Reference genome in fasta format"
            required = false
            arg_type = String
        "--reference-gff","-f"
            help = "Reference gff3 file"
            required = false
            arg_type = String
        "--use-precalculated-reference-transcripts-prefixes","-c"
            help = "Load polyA prefixes from previous run"
            action = :store_true
        "--incude-polyA-in-output","-i"
            help="Includes output polyA sequences as pseudo pair end's  into output fastq pair"
            action = :store_true
        "--processes","-p"
            help = "Number of additional julia workers for parallel procesing"
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
    eval(macroexpand(quote @everywhere using FMIndexes end))
    eval(macroexpand(quote @everywhere using IterTools end))
    eval(macroexpand(quote @everywhere push!(LOAD_PATH, ".") end))
    eval(macroexpand(quote @everywhere using PolyAAnalysis end))

    polyA_prefixes = get_polyA_prefixes_from_file(parsed_args["reference-transcripts"],
    parsed_args["reference-genome"],
    parsed_args["reference-gff"],
    minimum_not_polyA=parsed_args["minimum-length"],
    minimum_polyA_length=parsed_args["minimum-polyA-length"],
    number_of_workers=parsed_args["processes"],
    use_cached_results=parsed_args["use-precalculated-reference-transcripts-prefixes"])
    tic()

    trim_polyA_from_files(
        parsed_args["fastq-f"],
        parsed_args["fastq-r"],
        polyA_prefixes,
        parsed_args["processes"],
        parsed_args["output"],
        parsed_args["minimum-length"],
        parsed_args["minimum-polyA-length"],
        1, # maximum_non_A_symbols::Int64,
        3,  # maximum_distance_with_prefix_database
        4, # minimum_poly_A_between
        parsed_args["incude-polyA-in-output"]
        )
    toc()
end

main(ARGS)
