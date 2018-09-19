##
##  function detect_polyA_in_a_string(
##    fq_seq::String,
##    minimum_polyA_length::Int64,
##    maximum_non_A_symbols::Int64;
##    maximum_search_fragment_length::Int64=50,
##    debug=false
##    )::Bool
##

minimum_polyA_length=10
min_l = Int64(minimum_polyA_length / 2)
re = Regex("(A+[GTC])?(A{$min_l,})([GTC]A+)?")

seq = "AGAGAATTACAAATCAGAAGGGGAAGATTGACTGTTTAATAAAATGTGCTGGGAATCAATGGCAAATATACAATAAAATAAAATTATATTTTTATTCATATTATACCCCAAATAAATTTCAGATATACATAAAACCCAAAGTAAAATATT"
@test !detect_polyA_in_a_string(seq,10,1,re,debug=true)

seq = "GAATGTATGGTAGGAATGTATTCTCTTGTAGGAATGTAAATCTGTATTAAAAGGGGGTCCAAGCCAGGCCCCCAGGTCTTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAA"
@test detect_polyA_in_a_string(seq,10,1,re,debug=true)

seq = "AGAGAATTACAAATCAGAAGGGGAAGATTGACTGTTTAATAAAATGTGCTGGGAATCAATGGCAAATATACAATAAAATAAAATTATATTTTTATTCATATTATACCCCAAATAAATTTCAGATATACATAAAACCCAAAGTAAAATATT"
@test detect_polyA_in_a_string(seq,10,3,re,debug=true)

seq = "AGAGAATTACAAATCAGAAGGGGAAGATTGACTGTTTAATAAAATGTGCTGGGAATCAATGGCAAATATACAATAAAATAAAATTATATTTTTATTCATATTATACCCCAAATAAATTTCAGATATACATAAAACCCAAAGTAAAATATT"
@test !detect_polyA_in_a_string(seq,10,1,re,debug=true)

seq = "TCAAAAAAATAATATGGTAATAATAATAAAAGCAGTGCCCAGAGAACAAGGCTTGTTGGCTGTTCCACCCCAGGGGGCCCCTTGCACAGGCGGTGCCATCTCTGCCTCCCAAAGCTCTAAGAGCCACTGTCCCCCATCCCAAGAGA"
@test !detect_polyA_in_a_string(seq,10,1,re,debug=true)

seq = "TCAAAAAAATAATATGGTAATAATAATAAAAGCAGTGCCCAGAGAACAAGGCTTGTTGGCTGTTCCACCCCAGGGGGCCCCTTGCACAGGCGGTGCCATCTCTGCCTCCCAAAGCTCTAAGAGCCACTGTCCCCCATCCCAAGAGA"
@test !detect_polyA_in_a_string(seq,10,1,re,debug=true)

seq = "GGGGCGCGGCGGCGCGGGGGGGGAAAAAAAAACCCCCCCCCCCCCCCGGGGGGCCCC"
@test !detect_polyA_in_a_string(seq,10,1,re,debug=true)

seq = "GGGGCGCGGCGGCGCGGGGGGGGAAAATAAAAACCCCCCCCCCCCCCCGGGGGGCCCC"
@test detect_polyA_in_a_string(seq,10,1,re,debug=true)

seq = "AAGAAAAAAACAAACAACCCCATCAAAAAGTGGGCAAAGGATATAAATAGACAC"
@test detect_polyA_in_a_string(seq,10,1,re,debug=true)

seq = "TCAGGAGATCGAGACCATCCTGGCTAACATGGTGAAACCCCGTCTCTACTAAAAATAACAAAAAAAAA"
@test detect_polyA_in_a_string(seq,10,1,re,debug=true)

seq = "TGTGTGTGTGTGTGGTAAAAAAAAAAGT"
@test detect_polyA_in_a_string(seq,10,1,re,maximum_search_fragment_length=20)

seq = "TGTGTGTGTGTGTGGTAAAAAAAAAA"
@test detect_polyA_in_a_string(seq,10,1,re,maximum_search_fragment_length=20)

seq = "TGTGTGTGTGTGTGGTAAAAAAAAATAGGGG"
@test detect_polyA_in_a_string(seq,10,1,re,maximum_search_fragment_length=20)

seq = "TGTGTGTGTGTGTGGATAAAAAAAAA"
@test detect_polyA_in_a_string(seq,10,1,re,maximum_search_fragment_length=20)

seq = "TGTGTGTGTGTGTGGAATAAAAAAAAA"
@test detect_polyA_in_a_string(seq,10,1,re,maximum_search_fragment_length=20)

seq = "TGTGTGTGTGTGTGGAAAAAAAAATAA"
@test detect_polyA_in_a_string(seq,10,1,re,maximum_search_fragment_length=20)

seq = "TGTGTGTGTGTGTGGAATAAAAAAAAATAA"
@test detect_polyA_in_a_string(seq,10,1,re,maximum_search_fragment_length=20)

seq = "AAAAAAAAAAAAAAAAAAAATGGAACGCAGGGCAGGAACTCGTATTTGGGGGGAGATGGG"
@test detect_polyA_in_a_string(seq,10,1,re,debug=true)

seq = "TGTGTGTGTGTGTGGAAAAAAAAAT"
@test !detect_polyA_in_a_string(seq,10,1,re,maximum_search_fragment_length=20)

seq = "TGTGTGTGTGTGTGGTAAAAAAAAA"
@test !detect_polyA_in_a_string(seq,10,1,re,maximum_search_fragment_length=20)

seq = "TGTGTGTGTGTGTGGAAAAAAAAA"
@test !detect_polyA_in_a_string(seq,10,1,re,maximum_search_fragment_length=20)

seq = "AAAAAAAAAAAAAAATGTGTGTGTGTGTGGAAAAAAAAA"
@test !detect_polyA_in_a_string(seq,10,1,re,maximum_search_fragment_length=20)

seq = "AAAAAAAAAAAAAAA"
@test detect_polyA_in_a_string(seq,10,1,re,maximum_search_fragment_length=20)

seq = "AAAAAAAAA"
@test !detect_polyA_in_a_string(seq,10,1,re,maximum_search_fragment_length=20)
