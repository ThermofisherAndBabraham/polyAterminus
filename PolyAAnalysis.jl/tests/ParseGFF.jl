
# function get_transcripts_from_gff(record::BioSequences.FASTA.Record,
#    transDict::Dict)::Array{BioSequences.FASTA.Record,1}


dcAll = Dict("Seq" => [String("1-5,8-12,+"), String("1-5,+"), String("1-20,+")])
append!(dcAll["Seq"],[String("1-5,8-12,-"), String("1-5,-"), String("1-20,-")])
outRecordsAll=Array{BioSequences.FASTA.Record,1}()

# + strand testing
record = FASTA.Record("Seq", dna"GTGGTGGCGCTCACCGACTG")
dc = Dict("Seq" => [String("1-5,8-12,+")])
recordOut = FASTA.Record("Seq", dna"GTGGTCGCTC")
outRecords = Array{BioSequences.FASTA.Record,1}()
push!(outRecords,recordOut)
push!(outRecordsAll,recordOut)
@test outRecords == get_transcripts_from_dict(record,dc)

dc = Dict("Seq" => [String("1-5,+")])
recordOut = FASTA.Record("Seq", dna"GTGGT")
outRecords = Array{BioSequences.FASTA.Record,1}()
push!(outRecords,recordOut)
push!(outRecordsAll,recordOut)
@test outRecords == get_transcripts_from_dict(record,dc)

dc = Dict("Seq" => [String("1-20,+")])
recordOut = FASTA.Record("Seq", dna"GTGGTGGCGCTCACCGACTG")
outRecords = Array{BioSequences.FASTA.Record,1}()
push!(outRecords,recordOut)
push!(outRecordsAll,recordOut)
@test outRecords == get_transcripts_from_dict(record,dc)

# - strand testing
dc = Dict("Seq" => [String("1-5,8-12,-")])
recordOut = FASTA.Record("Seq", dna"GAGCGACCAC")
outRecords = Array{BioSequences.FASTA.Record,1}()
push!(outRecords,recordOut)
push!(outRecordsAll,recordOut)
@test outRecords == get_transcripts_from_dict(record,dc)

dc = Dict("Seq" => [String("1-5,-")])
recordOut = FASTA.Record("Seq", dna"ACCAC")
outRecords = Array{BioSequences.FASTA.Record,1}()
push!(outRecords,recordOut)
push!(outRecordsAll,recordOut)
@test outRecords == get_transcripts_from_dict(record,dc)

dc = Dict("Seq" => [String("1-20,-")])
recordOut = FASTA.Record("Seq", dna"CAGTCGGTGAGCGCCACCAC")
outRecords = Array{BioSequences.FASTA.Record,1}()
push!(outRecords,recordOut)
push!(outRecordsAll,recordOut)
@test outRecords == get_transcripts_from_dict(record,dc)

# Test all records in list
@test outRecordsAll == get_transcripts_from_dict(record,dcAll)
