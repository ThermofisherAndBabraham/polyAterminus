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
    trim_polyA_file_records,
    trim_polyA_from_fastq_record,
    get_polyA_prefixes,
    get_transcripts_from_gff,
    get_transcripts_from_dict,
    detect_polyA_in_a_string,
    extend_poly_A,
    trim_polyA_3end,
    check_polyA_prefixes,
    trim_polyA_from_fastq_pair,
    BamRead,
    Clust!,
    PolACalculus,
    ParseGFF3,
    ParseRecord,
    GetIntervalSet,
    ItsectCollection,
    rmdups,
    WrFrame

import BioAlignments: BAM
import DataFrames: DataFrame, DataFrameRow
import DataFrames: eachrow, deleterows!
import DataStructures
import CSV
import GenomicFeatures: GFF3
import GenomicFeatures: Interval
import GenomicFeatures: Strand
import GenomicFeatures: IntervalCollection
import GenomicFeatures: eachoverlap, isoverlapping, strand, metadata, seqname, first, last
import BioSequences: FASTQ
import BioSequences
import StringDistances: Levenshtein, evaluate
import TranscodingStreams
import CodecZlib
import BufferedStreams

include("MapPolyA.jl")
include("TrimmPolyA.jl")
include("ParseGFF.jl")

end  # module PolyAAnalysis
