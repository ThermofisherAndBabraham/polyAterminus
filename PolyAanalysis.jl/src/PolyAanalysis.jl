
# PolyAanalysis.jl
# =======
#
# A julia package for the representation and manipulation of biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

__precompile__()

module PolyAanalysis

export
    trim_polyA_file_records,
    trim_polyA_from_fastq_record,
    get_polyA_prefixes,
    get_polyA_prefixes_from_file

using BioSequences


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
    maximum_non_A_symbols::Int64=2,
    minimum_distance_from_non_poly_A::Int64=20)
    istream = FASTQ.Reader(open(input_fastq, "r"))
    ostream = FASTQ.Writer(open(output_fastq, "w"))

    ct_reads=0
    for fq in istream
        ct_reads+=1
        fqo,has_proper_polyA = trim_polyA_from_fastq_record(fq,
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
    minimum_distance_from_non_poly_A::Int64)::Tuple{FASTQ.Record,Bool}

    has_proper_polyA=false
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

    if has_proper_polyA
        return(fqo,has_proper_polyA)
    else
        return(fq,has_proper_polyA)
    end
end

"""
returns minimum kmer streches of not polyA ins a supplied transcript sequences
Arguments:
    sequence - string with polyA streches (type - BioSequences reference seq)
    minimum_not_polyA - minimum length of not polyA strech
    minimum_polyA_length - minimum length of polyA strech
    maximum_non_A_symbols - maximum numer of nonA symbols in polyA strech
    minimum_distance_from_non_poly_A - minimum number of any not polyA symbol in polyA strech from the non polyA
"""
function get_polyA_prefixes(sequence::BioSequences.ReferenceSequence,
    minimum_not_polyA::Int64,
    minimum_polyA_length::Int64,
    maximum_non_A_symbols::Int64,
    minimum_distance_from_non_poly_A::Int64)::Array{String,1}
    output=Array{String,1}()
    if maximum_non_A_symbols < 6 #beaks - not universal
        re=BioSequences.RE.Regex{DNA}("[ATGC]{$(minimum_not_polyA-1)}[TGC]A+[ATGC]A+[ATGC]A+[ATGC]A+[ATGC]A+")
    else
        error("Not posiible this value of errors")
    end

    for m in matchall(re, BioSequence{DNAAlphabet{4}}(sequence), false)
        polyA_part=m[minimum_not_polyA+1:end]
        polyA_length=length(polyA_part)

        #check if passes length requirement
        if polyA_length >= minimum_polyA_length
            #check if polyA mismatches are within limits

            # println((composition(polyA_part[1:minimum_polyA_length])[DNA_A] - minimum_polyA_length))
            if (composition(polyA_part[1:minimum_polyA_length])[DNA_A] - minimum_polyA_length) <= maximum_non_A_symbols
                # println(m, " all fragment")
                # println(polyA_part, " polyA_part")
                # println(polyA_length," ",minimum_polyA_length)
                #
                # println("-----------------")
                push!(output,m[1:minimum_not_polyA])
            end
        end

    end
    return(output)
end


"""
returns minimum kmer streches of not polyA ins a supplied transcript sequences from a file
Arguments:
    sequence - string with polyA streches (type - BioSequences reference seq)
    minimum_not_polyA - minimum length of not polyA strech
    minimum_polyA_length - minimum length of polyA strech
    maximum_non_A_symbols - maximum numer of nonA symbols in polyA strech
    minimum_distance_from_non_poly_A - minimum number of any not polyA symbol in polyA strech from the non polyA
"""
function get_polyA_prefixes_from_file(file::String,
    minimum_not_polyA::Int64,
    minimum_polyA_length::Int64,
    maximum_non_A_symbols::Int64,
    minimum_distance_from_non_poly_A::Int64)::Array{String,1}
    output=Array{String,1}()
    reader = FASTA.Reader(open(file, "r"))
    for record in reader
        rec_seq=ReferenceSequence(sequence(record))
        prefixes=get_polyA_prefixes(rec_seq,
        minimum_not_polyA,
        minimum_polyA_length,
        maximum_non_A_symbols,
        minimum_distance_from_non_poly_A)
        # if length(prefixes) > 0
        #     println(FASTA.identifier(record))
        #     println(prefixes)
        #     println(FASTA.sequence(record))
        # end
        output=vcat(output,prefixes)
    end
    close(reader)
    return(output)
end



end #module PolyAanalysis
