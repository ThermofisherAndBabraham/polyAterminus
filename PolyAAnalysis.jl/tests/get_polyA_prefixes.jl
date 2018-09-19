
##
## function get_polyA_prefixes(fasta_record::FASTA.Record,
##    minimum_not_polyA::Int64,
##    minimum_polyA_length::Int64)::Array{String,1}
##
minimum_not_polyA = 20
minimum_polyA_length = 10
re = Regex("([ATGC]{$minimum_not_polyA})A{$minimum_polyA_length,}")

# Test set 1 just a simple nt real example
seq = dna"ATCCTTGGCTTTCTTCCGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAA"
output = ["ATCCTTGGCTTTCTTCCGGG","CAAAAAAAAAAAAAAAAAAA", "GGGGGGGGGGGGGGGGGGGG"]
ref_seq=FASTA.Record("seq",seq)
@test get_polyA_prefixes(ref_seq,20,10,re)  == output

# Test set 2 no AAAA
seq = dna"ATCCTTGGCTTTCTTCCGGGCGGGGGGGGGGGGGGGGGGGGAAA"
output = []
ref_seq=FASTA.Record("seq",seq)
@test get_polyA_prefixes(ref_seq,20,10,re)  == output

# Test set 3
seq = dna"ATCCTTGGCTTTCTTCCGGGCGGGGGGGGGGGGGGGGAGAGAGAAAAAAAAAGAAAAAAAAAA"
output = ["GGGGAGAGAGAAAAAAAAAG"]
ref_seq=FASTA.Record("seq",seq)
@test get_polyA_prefixes(ref_seq,20,10,re)  == output

# Test set 4 Only AA...
seq = dna"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
output = ["AAAAAAAAAAAAAAAAAAAA"]
ref_seq=FASTA.Record("seq",seq)
@test get_polyA_prefixes(ref_seq,20,10,re)  == output

# Test set 5 Only AA... but too short
seq = dna"AAAAAAAAAAA"
output = []
ref_seq=FASTA.Record("seq",seq)
@test get_polyA_prefixes(ref_seq,20,10,re)  == output

# Test set 6
seq = dna"AAAAAAAAAAGGGGGGGGGGAAAAAAAAAA"
output = ["AAAAAAAAAAGGGGGGGGGG"]
ref_seq=FASTA.Record("seq",seq)
@test get_polyA_prefixes(ref_seq,20,10,re)  == output

# Test set 7 too short
seq = dna"GTCGTCGTC"
output = []
ref_seq=FASTA.Record("seq",seq)
@test get_polyA_prefixes(ref_seq,20,10,re)  == output
