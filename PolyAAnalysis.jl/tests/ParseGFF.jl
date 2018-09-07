
# function get_transcripts_from_gff(record::BioSequences.FASTA.Record,
#    transDict::Dict)::Array{BioSequences.FASTA.Record,1}

record = FASTA.Record("Seq", dna"GTGGTGGCGCTCACCGACTG")
dc = Dict("Seq" => [String("1-5,8-12+")])
recordOut = FASTA.Record("Seq", dna"GTGGTCGCTC")
outRecords = Array{BioSequences.FASTA.Record,1}()
push!(outRecords,recordOut)
@test outRecords == get_transcripts_from_gff!(record,dc)
