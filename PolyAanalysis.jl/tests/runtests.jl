module TestPolyAanalsis
using BioSequences
using Base.Test
push!(LOAD_PATH, "../../")
using PolyAanalysis

include("get_polyA_prefixes_from_file.jl")
exit()

@testset "PolyADatabaseconstractin" begin
include("get_polyA_prefixes.jl")
end


# @testset "PolyAtrimming" begin
# include("trim_polyA_from_fastq_record.jl")
# end

end #module TestPolyAanalsis
