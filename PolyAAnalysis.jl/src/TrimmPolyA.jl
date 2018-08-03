#!/usr/bin/env julia



"""
checks if a string of a read contains polyA strechA strating from the 3' end
    fq_seq - string for testing
    minimum_polyA_length - minimum length of polyA strech
    maximum_non_A_symbols - maximum numer of nonA symbols in polyA strech
"""
function detect_polyA_in_a_string(
    fq_seq::String,
    minimum_polyA_length::Int64,
    maximum_non_A_symbols::Int64;
    debug=false
    )::Bool
    has_poly_a=false
    fq_length=length(fq_seq)
    window_position_from_3_end=0
    for startini in minimum_polyA_length:fq_length
        window_position_from_3_end+=1
        start=fq_length-startini+1
        substring=String(fq_seq[start:start+minimum_polyA_length-1])
        if substring[1]=='A' #first symbol of apolyA strech must be A
            ct=0
            for symb in substring
                if symb != 'A'
                    ct+=1
                end
            end
            if ct <=maximum_non_A_symbols
                has_poly_a=true
                if debug
                    println("substring having less than $maximum_non_A_symbols non A symbols")
                    println(substring)
                end
                break
            end
        end
    end
    return(has_poly_a)
end

"""
Tries to merge and extend forward read
    fq_seq1 - string of forward read
    fq_seq2 - string of reverse read reverse complement
"""
function extend_poly_A(
    fq_seq1::FASTQ.Record,
    fq_seq2::FASTQ.Record)::FASTQ.Record
    return(fq_seq1)
end

"""
check if the read is in pefixes list of naturall polyA sreches
    fqo_trimmed - fastq record
    prefixes - array of natural prefixes
    maximum_distance_with_prefix_database - maximum Levenstain distance
"""

function check_polyA_prefixes(
    fqo_trimmed::FASTQ.Record,
    prefixes::Array{String,1},
    maximum_distance_with_prefix_database::Int64



    )::Bool
    has_no_match=true
    for pref in prefixes
        if  evaluate(Levenshtein(), String(FASTQ.sequence(fqo_trimmed)),pref) <= maximum_distance_with_prefix_database
            has_no_match=false
            break
        end
    end
    return(has_no_match)
end

"""
function that trimmes 3' end of polyA having sequence from non A symbols that might originate fue to reamins of adapters or sequencing artefacts
Arguments:
    seq - read sequence in String
    minimum_poly_A_between - minimum length of polyA strech that might occure between non A symbols forming a fragment to be trimmed off (starting from the very 3' end)
"""
function trim_polyA_3end(seq::String,
    minimum_poly_A_between::Int64
    )::Tuple{String,Int64}
    length_sequence=length(seq)
    i=length(seq)
    trimming_position=0
    polyA_strech_length=0
    while (polyA_strech_length <= minimum_poly_A_between) && (i>0)
        symbol=seq[i]
        if symbol != 'A'
            trimming_position=i
            polyA_strech_length=0
        else
            polyA_strech_length+=1
        end
        i-=1
    end
    if trimming_position==0
        out=(seq,0)
    else
        out=(seq[1:trimming_position-1],length_sequence-trimming_position+1)
    end
    return(out)
end

"""
counts not A in a string
"""
function count_nona(seq::String)::Int64
    ct=0
    for s in seq
        if s != 'A'
            ct+=1
        end
    end
    return(ct)
end

"""
finds first not A from the end (3')
"""
function first_nona(seq::String)::Int64
    ct=0
    l=length(seq)
    for i in 1:l
        pos=l-ct
        if seq[pos]!='A'
            break
        end
        ct+=1
    end
    return(ct)
end




"""
function that trimmes 3 polyA tail of a given fastq entry and outputs trimmed read with number of polyA indicated in read's name
Arguments:
    minimum_not_polyA - minimum length of not polyA strech in a read
    minimum_polyA_length - minimum length of polyA strech
    minimum_poly_A_between - minimum length of polyA strech that might occure between non A symbols forming a fragment to be trimmed off (starting from the very 3' end)
Keywords arguments:
    debug - turns out additional output
    window_length - windows size for polyA trimming
    max_nonA_in_window - tollerated nonA number per window

The algorithms works like:
    1. Trimmes offs remains of adapters at the 3'
    2. Search for a window wich exceeds nonA limit and has no other high A content window towards 5' end
    3. Ajusts the A position in the window
"""
function trim_polyA_from_fastq_record(fq::FASTQ.Record,
    minimum_not_polyA::Int64,
    minimum_polyA_length::Int64,
    minimum_poly_A_between::Int64;
    debug=false,
    window_length=10,
    max_nonA_in_window=1,
    limit_for_3end_trimming=20
    )::Tuple{FASTQ.Record,Bool,Bool}
    has_proper_polyA=false
    polyA_detected=true
    seq_ini=String(FASTQ.sequence(fq),)
    seq,trimmed_to_polyA=trim_polyA_3end(seq_ini,minimum_poly_A_between)

    seq_len=length(seq)
    name=FASTQ.identifier(fq)
    description=FASTQ.description(fq)
    quality=FASTQ.quality(fq)[1:seq_len]
    i=seq_len
    polyA_start=0
    first_nona_position=0
    if debug
        println("-----------------------------------------------------------------------------------------------------------------------------\nSequence for trimming:\n$seq_ini")
        println("Sequence with chopped 3'end adapter remains:\n$seq")
        if !(polyA_detected)
            println("WARNING! this sequence will be concidered as not having polyA tail due to trimming of 3end adapter remains too much (> $limit_for_3end_trimming bp)")
        end
    end

    # search for a polyA start with a sliding window from the 3' end
    #Analysis limit to the cases where we meet a new non A symbol (going from 3' end)
    initial_position=0 # position of polyA start based on windows.
    first_nona_position=0
    while i > 0 & initial_position ==0
        s=seq[i]
        if s != 'A'
            fragment_start=i
            fragment_end=i+window_length-1
            # println("fragment_start:fragment_end-ini seq_len $fragment_start $fragment_end $seq_len")
            # println(seq)
            if fragment_end > seq_len #handle the cases when there are non A at the 3' end
                fragment_end_ajusted=seq_len
                fragment_start_ajusted=i - (fragment_end-seq_len)
                fragment_start=fragment_start_ajusted
                fragment_end=fragment_end_ajusted
            end

            if fragment_start < 1
                fragment_start=1
            end
            # println("fragment_start:fragment_end $fragment_start $fragment_end")
            fragment=seq[fragment_start:fragment_end]
            number_of_nonA=count_nona(fragment)
            if debug
                println("Fragment with non A, starting position $i, with $number_of_nonA nonA symbols")
                println(fragment)
            end
            if max_nonA_in_window < number_of_nonA
                #test if a window towards 5' also contains non A sequences
                to_5end_frag_end=fragment_start-1
                to_5end_frag_start=to_5end_frag_end-window_length+1
                if to_5end_frag_start < 1
                     to_5end_frag_start = 1
                end
                if to_5end_frag_start == 1 || (count_nona(seq[to_5end_frag_start:to_5end_frag_end]) > max_nonA_in_window)
                    initial_position=fragment_end
                    first_nona_position=initial_position-first_nona(seq[1:initial_position])
                    if debug
                        println("Initial non A sequence:")
                        println(seq[1:initial_position])
                        println("Adjusted non A sequence:")
                        println(seq[1:first_nona_position])
                    end
                    break
                end
            end
        end
        i-=1
    end

    polyA_start=first_nona_position+1
    polyA_length=seq_len-polyA_start+1

    if (polyA_start > 0) & (polyA_start > minimum_not_polyA)
        not_poly_a_seq=seq[1:polyA_start-1]
        nname=string(polyA_length)*":A:"*name
        nquality=quality[1:polyA_start-1]
        #fqo=FastqRead(nname, dna(not_poly_a_seq), strand, nquality)
        fqo=FASTQ.Record(nname, description, not_poly_a_seq, nquality)
        has_proper_polyA=true
    end

    if (trimmed_to_polyA > limit_for_3end_trimming) || (polyA_length < minimum_polyA_length)
        polyA_detected=false
    end

    if has_proper_polyA
        return(fqo,has_proper_polyA,polyA_detected)
    else
        return(fq,has_proper_polyA,polyA_detected)
    end

end


"""
returns minimum kmer streches of not polyA ins a supplied transcript sequences
Arguments:
    sequence - string with polyA streches (type - BioSequences reference seq)
    minimum_not_polyA - minimum length of not polyA strech
    minimum_polyA_length - minimum length of polyA strech
"""
function get_polyA_prefixes(fasta_record::BioSequences.FASTA.Record,
    minimum_not_polyA::Int64,
    minimum_polyA_length::Int64)::Array{String,1}
    output=Array{String,1}()
    re=Regex("[ATGC]{$minimum_not_polyA}A{$minimum_polyA_length}")
    for m in eachmatch(re, String(BioSequences.sequence(fasta_record)))
        prefix_part=m.match[1:minimum_not_polyA]
        push!(output,prefix_part)
    end
    return(output)
end

"""
returns minimum kmer streches of not polyA ins a supplied transcript sequences
Arguments:
    sequence - string with polyA streches (type - BioSequences reference seq)
    minimum_not_polyA - minimum length of not polyA strech
    minimum_polyA_length - minimum length of polyA strech
    counter - shared array for tracking processed reads - length matches julia workers
"""
function get_polyA_prefixes(fasta_record::BioSequences.FASTA.Record,
    minimum_not_polyA::Int64,
    minimum_polyA_length::Int64,
    counter::SharedArray)::Array{String,1}
    output=get_polyA_prefixes(fasta_record,minimum_not_polyA,minimum_polyA_length)
    counter[(myid()-1)] += 1
    total_jobs = sum(counter)
    if mod(total_jobs,10000)==0
      print(STDERR,"Number of parsed transcripts = $total_jobs\r")
    end
    return(output)
end

"""
finds and trims polyA having reads from a pair of fastq records
Arguments:
    fastq1 - FASTQ record
    fastq2 - FASTQ record
    prefixes - array of prefixes of natural polyA
    minimum_not_polyA - minimum length of not polyA strech in a read
    minimum_polyA_length - minimum length of polyA strech
    maximum_non_A_symbols - maximum numer of nonA symbols in polyA strech
    maximum_distance_with_prefix_database - allowed Levenshtein distance between prefix of natural polyA in transcripts and the read
    minimum_poly_A_between - minimum length of polyA strech that might occure between non A symbols forming a fragment to be trimmed off (starting from the very 3' end)

"""

function trim_polyA_from_fastq_pair(
    fastq1::FASTQ.Record,
    fastq2::FASTQ.Record,
    prefixes::Array{String,1},
    minimum_not_polyA::Int64,
    minimum_polyA_length::Int64,
    maximum_non_A_symbols::Int64,
    maximum_distance_with_prefix_database::Int64,
    minimum_poly_A_between::Int64;
    debug_id="None",
    debug=false
    )
    polyA_detected=false
    has_proper_polyA=false
    revseq_rev_read=String(BioSequences.reverse_complement!(FASTQ.sequence(fastq2)))
    seq_for_read=String(FASTQ.sequence(fastq1))

    if debug
        println("seq_for_read", seq_for_read)
        println("revseq_rev_read", revseq_rev_read)
        println(detect_polyA_in_a_string(seq_for_read,minimum_polyA_length,maximum_non_A_symbols))
        println(detect_polyA_in_a_string(revseq_rev_read,minimum_polyA_length,maximum_non_A_symbols))
    end


    if detect_polyA_in_a_string(seq_for_read,minimum_polyA_length,maximum_non_A_symbols) &&
        detect_polyA_in_a_string(revseq_rev_read,minimum_polyA_length,maximum_non_A_symbols)
            #tries sto esxtend polyA
            cosensus_extended_fwd_read=extend_poly_A(fastq1,fastq2)
            fqo_trimmed, has_proper_polyA, polyA_detected = trim_polyA_from_fastq_record(cosensus_extended_fwd_read,
                minimum_not_polyA,
                minimum_polyA_length,
                minimum_poly_A_between,debug=debug)
            if debug
                println("Output of $trim_polyA_from_fastq_record")
                println("Input seq",cosensus_extended_fwd_read)
                println(fqo_trimmed, has_proper_polyA, polyA_detected)
                println(minimum_not_polyA)
                println(minimum_polyA_length)
                println(minimum_poly_A_between)
            end

            #check if the read is in pefixes list of natural polyA streches
            if has_proper_polyA
                has_proper_polyA = check_polyA_prefixes(fqo_trimmed,prefixes,maximum_distance_with_prefix_database)
                if !has_proper_polyA
                    fqo_trimmed=cosensus_extended_fwd_read
                end
            end
    else
        fqo_trimmed=fastq1
    end
    return(fqo_trimmed, has_proper_polyA, polyA_detected)
end



"""
finds and trims polyA having reads from a pair of fastq records
Arguments:
    fastq1 - FASTQ record
    fastq2 - FASTQ record
    prefixes - array of prefixes of natural polyA
    minimum_not_polyA - minimum length of not polyA strech in a read
    minimum_polyA_length - minimum length of polyA strech
    maximum_non_A_symbols - maximum numer of nonA symbols in polyA strech
    maximum_distance_with_prefix_database - allowed Levenshtein distance between prefix of natural polyA in transcripts and the read
    minimum_poly_A_between - minimum length of polyA strech that might occure between non A symbols forming a fragment to be trimmed off (starting from the very 3' end)
    include_polyA - Includes output polyA sequences as pseudo pair end's  into output fastq
    fastqo1 - output stream for forward fastq.gz
    fastqo2 - output stream for reverse fastq.gz
    fastqo_s - output stream  for polyA trimmed and marked reads (forward or merged reads only)
    fastqo_d - output stream  for discarded reads (forward or merged reads only)
    debug_id - name of fastq entry that should be debuged and additional infor printed our
"""

function trim_polyA_from_fastq_pair(
    fastq1::FASTQ.Record,
    fastq2::FASTQ.Record,
    prefixes::Array{String,1},
    minimum_not_polyA::Int64,
    minimum_polyA_length::Int64,
    maximum_non_A_symbols::Int64,
    maximum_distance_with_prefix_database::Int64,
    minimum_poly_A_between::Int64,
    include_polyA::Bool,
	ct_all::SharedArray,
	ct_pair_without_polyA::SharedArray,
	ct_pair_with_proper_polyA::SharedArray,
	ct_pair_with_discarded_polyA::SharedArray;
    debug_id="None",
    debug=false
    )
    #initial filtering
    if debug_id == FASTQ.identifier(fastq1)
        debug=true
    else
        debug=false
    end

	ct_all[(myid()-1)] += 1


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
    #detect_polyA_in_fastq_record(records[1],minimum_polyA_length,maximum_non_A_symbols)

	#get global variables of output streams
	global fastqo1 = Main.fastqo1
	global fastqo2 = Main.fastqo2
	global fastqo_d = Main.fastqo_d
	global fastqo_s = Main.fastqo_s

    if !polyA_detected
		ct_pair_without_polyA[(myid()-1)] += 1
		#write(test,fastq1)
		println(fastqo1,fastq1)
        println(fastqo2,fastq2)
		flush(fastqo1)
		flush(fastqo2)
    else
        if has_proper_polyA
            println(fastqo_s,fqo_trimmed)
			flush(fastqo_s)
			ct_pair_with_proper_polyA[(myid()-1)] += 1
            #get reverse complement of the polyA read
            if include_polyA
                name=FASTQ.identifier(fqo_trimmed)
                description=FASTQ.description(fqo_trimmed)
                rev_quality=reverse(FASTQ.quality(fqo_trimmed))
                rev_seq=BioSequences.reverse_complement!(FASTQ.sequence(fqo_trimmed))
                fqo_trimmed_rev=FASTQ.Record(name, description, rev_seq, rev_quality)
				println(fastqo1, fqo_trimmed)
                println(fastqo2, fqo_trimmed_rev)
				flush(fastqo1)
				flush(fastqo2)
            end
        else

            println(fastqo_d,fqo_trimmed)
			flush(fastqo_d)
			ct_pair_with_discarded_polyA[(myid()-1)] += 1
        end
    end
	if mod(sum(ct_all),1000)==0
		println("Parsed reads: ",sum(ct_all),
		" without polyA: ",sum(ct_pair_without_polyA),
		" with proper polyA: ",sum(ct_pair_with_proper_polyA),
		" with discarded polyA: ",sum(ct_pair_with_discarded_polyA),"\r")
	end

end
