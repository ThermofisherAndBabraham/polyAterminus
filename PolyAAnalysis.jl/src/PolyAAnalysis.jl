#!/usr/bin/env julia

# PolyAAnalysis.jl
# =======
#
# A julia package for the representation and manipulation of biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

__precompile__()

module PolyAAnalysis

export
    annotate_polya_sites,
    annotate_clusters,
    check_polyA_prefixes,
    clust!,
    cluster!,
    detect_polyA_in_a_string,
    enumeratenames!,
    extend_poly_A,
    get_polyA_prefixes,
    get_split_key,
    get_transcripts_from_dict,
    get_transcripts_from_gff,
    getintervals,
    parseGFF3,
    parserecord,
    readbam,
    rmdups,
    stats_poly_a,
    stats_clusters,
    trim_polyA_3end,
    trim_polyA_file_records,
    trim_polyA_from_fastq_pair,
    trim_polyA_from_fastq_record

import BioAlignments: BAM
import BioSequences: BioSymbols, @dna_str, FASTA, reverse_complement!, sequence
import BioSequences: FASTQ
import BufferedStreams
import CodecZlib
import CSV
import DataFrames: DataFrame, DataFrameRow
import DataFrames: eachrow, deleterows!
import DataStructures
import GenomicFeatures: eachoverlap, isoverlapping, strand, metadata, seqname, first, last
import GenomicFeatures: GFF3
import GenomicFeatures: Interval
import GenomicFeatures: IntervalCollection
import GenomicFeatures: Strand
import StringDistances: Levenshtein, evaluate
import TranscodingStreams
import FMIndexes: FMIndexes
import StatsBase: nquantile


include("AnnotatePolyA.jl")
include("ParseGFF.jl")
include("TrimmPolyA.jl")

end  # module PolyAAnalysis
