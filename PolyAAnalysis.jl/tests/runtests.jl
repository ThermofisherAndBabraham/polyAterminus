module TestPolyAAnalsis

using BioSequences
using Base.Test
push!(LOAD_PATH, "../../")
using DataFrames
using DataStructures
import GenomicFeatures: GFF3
import GenomicFeatures: Interval
import GenomicFeatures: Strand
import GenomicFeatures: IntervalCollection
using PolyAAnalysis

@testset "PolyAtrimming" begin
include("trim_polyA_3end.jl")
include("trim_polyA_from_fastq_record.jl")
include("detect_polyA_in_a_string.jl")
include("get_polyA_prefixes.jl")
end

@testset "AnnotatePolyA" begin
include("AnnotatePolyA.jl")
end

@testset "ParseGFF" begin
include("ParseGFF.jl")
end
end #module TestPolyAanalsis
