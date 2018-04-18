module TestPolyAanalsis
using BioSequences
using Base.Test
push!(LOAD_PATH, "../../")
using PolyAanalysis

@testset "PolyAtrimming" begin
include("trim_polyA_from_fastq_record.jl")
end


@testset "PolyADatabase" begin
include("get_polyA_prefixes.jl")
include("get_polyA_prefixes_from_file.jl")
end





end #module TestPolyAanalsis
