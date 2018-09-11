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
    BamRead,
    check_polyA_prefixes,
    Clust!,
    detect_polyA_in_a_string,
    extend_poly_A,
    get_polyA_prefixes,
    get_transcripts_from_dict,
    get_transcripts_from_gff,
    GetIntervalSet,
    ItsectCollection,
    ParseGFF3,
    ParseRecord,
    PolACalculus,
    rmdups,
    trim_polyA_3end,
    trim_polyA_file_records,
    trim_polyA_from_fastq_pair,
    trim_polyA_from_fastq_record,
    WrFrame

import BioAlignments: BAM
import BioSequences
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

include("MapPolyA.jl")
include("TrimmPolyA.jl")
include("ParseGFF.jl")

end  # module PolyAAnalysis
