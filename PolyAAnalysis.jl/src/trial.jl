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
        println(start," ",start+minimum_polyA_length-1, " ",fq_length)
        substring=String(fq_seq[start:start+minimum_polyA_length-1])
        println(substring)
        println(substring)
        if substring[1]=='A' #first symbol must be A
            ct=0
            for symb in substring
                if symb != 'A'
                    ct+=1
                end
            end

            if ct <=maximum_non_A_symbols
                has_poly_a=true
                println(ct,substring)
                break
            end
        end
    end
    return(has_poly_a)
end

detect_polyA_in_a_string("ATAAAAAAAAAGGGAGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGAAGAAAAAAAAAAAAAAAAAAGGGGGGGGGC",20,3)
