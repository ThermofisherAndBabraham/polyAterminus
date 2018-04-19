module TestPolyAanalsis

using BioSequences
using Base.Test
push!(LOAD_PATH, "../../")
using PolyAAnalysis


@testset "PolyAtrimming" begin
include("trim_polyA_from_fastq_record.jl")
end

@testset "MapPolyA" begin
include("MapPolyA_tests.jl")
end

end #module TestPolyAanalsis
