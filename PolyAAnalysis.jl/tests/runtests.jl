module TestPolyAAnalsis

using BioSequences
using Base.Test
push!(LOAD_PATH, "../../")
using DataFrames
using PolyAAnalysis

@testset "PolyAtrimming" begin
include("trim_polyA_from_fastq_record.jl")
include("trim_polyA_from_fastq_record.jl")
end

@testset "MapPolyA" begin
include("MapPolyA.jl")
end



end #module TestPolyAanalsis
