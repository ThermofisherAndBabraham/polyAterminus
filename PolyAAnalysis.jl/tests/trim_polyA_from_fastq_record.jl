# function trim_polyA_from_fastq_record(fq::FASTQ.Record,
#     minimum_not_polyA::Int64,
#     minimum_polyA_length::Int64,
#     maximum_non_A_symbols::Int64,
#     minimum_distance_from_non_poly_A::Int64)::Tuple{FASTQ.Record,Bool}


#test set 1
seq = "ATCCTTGGCTTTCTTCCGGGCCGGTTGAAAAAAAAAAAAAAAAAAAAAAAAAAA"
trimmed = "ATCCTTGGCTTTCTTCCGGGCCGGTTG"
description="blabla"
len = 27
q=fill(20,length(seq))
q_trimmed=fill(20,length(trimmed))
fq=FASTQ.Record("testing_seq",description,seq,q)
fq_trimmed=FASTQ.Record("$len:A:testing_seq",description,trimmed,q_trimmed)
a=trim_polyA_from_fastq_record(fq,20,10,1,1)
#println()

#test trimming
@test FASTQ.sequence(trim_polyA_from_fastq_record(fq,20,10,1,1)[1])  == FASTQ.sequence(fq_trimmed)
#test length calculation
@test FASTQ.identifier(trim_polyA_from_fastq_record(fq,20,10,1,1)[1])  == FASTQ.identifier(fq_trimmed)
