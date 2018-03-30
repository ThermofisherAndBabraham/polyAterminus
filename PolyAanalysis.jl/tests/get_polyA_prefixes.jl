# """
# returns minimum kmer streches of not polyA ins a supplied transcript sequences
# Arguments:
#     sequence - string with polyA streches
#     minimum_not_polyA - minimum length of not polyA strech
#     minimum_polyA_length - minimum length of polyA strech
#     maximum_non_A_symbols - maximum numer of nonA symbols in polyA strech
#     minimum_distance_from_non_poly_A - minimum number of any not polyA symbol in polyA strech from the non polyA
# """
# function get_polyA_prefixes(sequence::String,minimum_not_polyA::Int64,minimum_polyA_length::Int64,maximum_non_A_symbols::Int64)::Array{String,1}
#




#test set 1 just a simple nt real example
seq = "ATCCTTGGCTTTCTTCCGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAA"
output = ["ATCCTTGGCTTTCTTCCGGG", "GGGGGGGGGGGGGGGGGGGG"]
ref_seq=ReferenceSequence(BioSequence{DNAAlphabet{4}}(seq))

#println()

#test trimming
@test get_polyA_prefixes(ref_seq,20,10,1,10)  == output
