#!/usr/bin/env julia

#= Custom library =#
using ArgParse.ArgParseSettings    #=                                    =#
using ArgParse.@add_arg_table      #= For parsing command-line options   =#
using ArgParse.parse_args          #=                                    =#
using OpenGene




function trim(output_fastq::String,
    input_fastq::String,
    minimum_not_polyA::Int64;
    minimum_polyA_length::Int64=50,
    maximum_non_A_symbols::Int64=2,
    minimum_distance_from_non_poly_A::Int64=20)

    istream = fastq_open(input_fastq)

    ostream = fastq_open(output_fastq,"w")

    testname="ST-E00243:413:HKCCHCCXY:6:2209:27184:47070" #for testing read name
    ct_reads=0
    while (fq = fastq_read(istream))!=false
        ct_reads+=1
        seq=fq.sequence.seq
        name=fq.name
        quality=fq.quality
        strand=fq.strand
        seq_len=length(seq)
        i=seq_len
        ct_all_symbols=0
        ct_A=0
        ct_not_A=0
        last_not_ca_position=seq_len
        polyA_start=0 #search for it in the loop #iterate from the 3' end
        # if contains(name,testname)
        #     println(seq)
        # end

        while (ct_not_A < maximum_non_A_symbols) & (i > 0  )
            s=seq[i]::Char
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
            if contains(name,testname)
                println("i s polyA_start ct_not_A $i $s $polyA_start $ct_not_A ")
            end
            i-=1
        end


        polyA_length=seq_len-polyA_start+1
        if (polyA_start > 0) & (polyA_start > minimum_not_polyA)
            not_poly_a_seq=seq[1:polyA_start-1]
            nname=name[1:1]*string(polyA_length)*"_"*name[2:end]
            nquality=quality[1:polyA_start-1]
            #fqo=FastqRead(nname, dna(not_poly_a_seq), strand, nquality)
            fqo=FastqRead(nname, dna(seq), strand, quality) #for testing
            fastq_write(ostream, fqo)
            if contains(name,testname)
                println("-----------------------------------")
                println(polyA_start > 0)
                println(seq_len, " ", polyA_start, " ", minimum_not_polyA)
                println(nname," ")
                println(not_poly_a_seq)
                println(seq)
            end
        end


        if mod(ct_reads,1000)==0
          print(STDERR, "PROCESSED $ct_reads FASTQ records \r")
        end

    end

    close(ostream)

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
    #polymorphism_limit::Float64=0.1, minimum_coverage::Int64=5 , mutations_number_limit
     end

    parsed_args = parse_args(arg_parse_settings) #= In order to use variables, =#
                                                 #= use parsed_args["foo_bar"] =#
    #= Main code =#
    trim(parsed_args["output"],
    parsed_args["input"],
    parsed_args["minimum-length"])
end

main(ARGS)
