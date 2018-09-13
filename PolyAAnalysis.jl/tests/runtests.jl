module TestPolyAAnalsis

using BioSequences
using Base.Test
push!(LOAD_PATH, "../../")
using DataFrames
import GenomicFeatures: GFF3
import GenomicFeatures: Interval
import GenomicFeatures: Strand
import GenomicFeatures: IntervalCollection
using PolyAAnalysis

@testset "PolyAtrimming" begin
include("trim_polyA_3end.jl")
include("trim_polyA_from_fastq_record.jl")
end

@testset "AnnotatePolyA" begin
include("AnnotatePolyA.jl")
end

@testset "ParseGFF" begin
include("ParseGFF.jl")
end
end #module TestPolyAanalsis
