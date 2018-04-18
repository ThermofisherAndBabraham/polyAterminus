using BioSequences
using Base.Test
push!(LOAD_PATH, "../../")
using PolyAanalysis
# """
# returns minimum kmer streches of not polyA ins a supplied transcript sequences
# Arguments:
#     file - fastafile with transcripts
#     minimum_not_polyA - minimum length of not polyA strech
#     minimum_polyA_length - minimum length of polyA strech
#     maximum_non_A_symbols - maximum numer of nonA symbols in polyA strech
#     minimum_distance_from_non_poly_A - minimum number of any not polyA symbol in polyA strech from the non polyA
# """
# function get_polyA_prefixes(sequence::String,minimum_not_polyA::Int64,minimum_polyA_length::Int64,maximum_non_A_symbols::Int64)::Array{String,1}
#

test_file="testing_datasets/transcripts.fa"
#test trimming
@test String["ATCCTTGGCTTTCTTCCGGG","GGCGGGGAAGGGGGGGGGGG"] == get_polyA_prefixes_from_file(test_file,20,40,1,20)
