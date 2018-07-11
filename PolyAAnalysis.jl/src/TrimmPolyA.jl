#!/usr/bin/env julia


"""
Trims a R1 fastq file and writes polyA reads to output file

Arguments:
    output_fastq - output fastq file name
    input_fastq - input fastq file name
    minimum_not_polyA - minimum length of not polyA strech in a read
    minimum_polyA_length - minimum length of polyA strech
    maximum_non_A_symbols - maximum numer of nonA symbols in polyA strech
    minimum_distance_from_non_poly - minimum distance of non A symbols in polyA strech from the nonA sequence
"""
function trim_polyA_file_records(output_fastq::String,
    input_fastq::String,
    minimum_not_polyA::Int64;
    minimum_polyA_length::Int64=50,
    maximum_non_A_symbols::Int64=1,
    minimum_distance_from_non_poly_A::Int64=20)
    istream = BioSequences.FASTQ.Reader(open(input_fastq, "r"))
    ostream = BioSequences.FASTQ.Writer(open(output_fastq, "w"))
    #prepare polyA streches database from transcripts

    ct_reads=0
    for fq in istream
        ct_reads+=1
        fqo,has_proper_polyA, polyA_detected = trim_polyA_from_fastq_record(fq,
            minimum_not_polyA,
            minimum_polyA_length,
            maximum_non_A_symbols,
            minimum_distance_from_non_poly_A)
        if has_proper_polyA
            write(ostream,fqo)
        end
        if mod(ct_reads,1000)==0
          print(STDERR, "PROCESSED $ct_reads FASTQ records \r")
        end

    end

    close(istream)
    close(ostream)

end



"""
checks if a string of a read contains polyA strechA strating from the 3' end
    fq_seq - string for testing
    minimum_polyA_length - minimum length of polyA strech
    maximum_non_A_symbols - maximum numer of nonA symbols in polyA strech
    maximum_number_of_adapter_remains_3_end - some adapter remains can be left after trimming
"""
function detect_polyA_in_a_string(
    fq_seq::String,
    minimum_polyA_length::Int64,
    maximum_non_A_symbols::Int64;
    maximum_number_of_adapter_remains_3_end=10)::Bool
    has_poly_a=false
    fq_length=length(fq_seq)
    window_position_from_3_end=0
    for startini in minimum_polyA_length:fq_length
        window_position_from_3_end+=1
        if window_position_from_3_end > maximum_number_of_adapter_remains_3_end
            break
            #seach of polyA strech discontinued after sertain length after 3' end
            #some adapter remains can be left after trimming
        end
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
function that trimmes 3 polyA tail of a given fastq entry and outputs trimmed read with number of polyA indicated in read's name
Arguments:
    minimum_not_polyA - minimum length of not polyA strech in a read
    minimum_polyA_length - minimum length of polyA strech
    maximum_non_A_symbols - maximum numer of nonA symbols in polyA strech
    minimum_distance_from_non_poly - minimum distance of non A symbols in polyA strech from the nonA sequence
"""
function trim_polyA_from_fastq_record(fq::FASTQ.Record,
    minimum_not_polyA::Int64,
    minimum_polyA_length::Int64,
    maximum_non_A_symbols::Int64,
    minimum_distance_from_non_poly_A::Int64)::Tuple{FASTQ.Record,Bool,Bool}

    has_proper_polyA=false
    polyA_detected=false
    seq=String(FASTQ.sequence(fq))
    name=FASTQ.identifier(fq)
    description=FASTQ.description(fq)
    quality=FASTQ.quality(fq)
    seq_len=length(seq)
    i=seq_len
    ct_all_symbols=0
    ct_A=0
    ct_not_A=0
    last_not_ca_position=seq_len
    polyA_start=0

    while (ct_not_A < maximum_non_A_symbols) & (i > 0  )
        s=seq[i]
        ct_all_symbols+=1
        curr_symbol_is_A::Bool=false
        if s=='A'
            ct_A+=1
            curr_symbol_is_A=true
        else
            ct_not_A+=1
            last_not_ca_position=i
        end

        if curr_symbol_is_A && ((seq_len-i) >= minimum_polyA_length) && ((last_not_ca_position-i) > minimum_distance_from_non_poly_A)
            polyA_start=i
        end
        i-=1
    end

    polyA_length=seq_len-polyA_start+1
    if (polyA_start > 0) & (polyA_start > minimum_not_polyA)
        not_poly_a_seq=seq[1:polyA_start-1]
        nname=string(polyA_length)*":A:"*name
        nquality=quality[1:polyA_start-1]
        #fqo=FastqRead(nname, dna(not_poly_a_seq), strand, nquality)
        fqo=FASTQ.Record(nname, description, not_poly_a_seq, nquality)
        has_proper_polyA=true
    end

    polyA_detected=(polyA_start > 0)

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
