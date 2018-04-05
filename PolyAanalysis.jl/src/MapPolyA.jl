#!/usr/bin/env julia

using ArgParse.ArgParseSettings
using ArgParse.@add_arg_table
using ArgParse.parse_args
using Bio.Align
using DataFrames
using DataStructures
using CSV
using Bio.Intervals

"""
"Turi buti: PolyA parsing is read name -> clustering -> statistics -> gtf, arba bed."
"Count total reads, aligned reads in intervals and general coverage."
"""
function BamRead(bam::String)

    reader = open(Bio.Align.BAMReader, bam)
    record = Bio.Align.BAMRecord()
    treads = Int64(0)
    passreads = Int64(0)
    pareads = Int64(0)
    pclus = Dict{String,Array{Int16,1}}()

    while !eof(reader)
        read!(reader, record)
        rfname::String = refname(record)
        flag = Bio.Align.flag(record)
        treads += 1

        if (flag&4 == 0) && (flag&256 == 0) && (flag&2048 == 0)  # filter out only mapped and primary alignments
            spl = split(seqname(record),":")

            if spl[2] == "A"
                pareads += 1
                rpos = rightposition(record)
                chrpos = rfname * "::" * string(rpos)
                pclus = PolAClus(pclus, chrpos, parse(Int16, spl[1]))
            end
            passreads += 1
        end
    end

    return treads, passreads, pareads, pclus
end


function PolAClus(d::Dict{String,Array{Int16,1}}, p::String, l::Int16)

    if haskey(d, p)
        v = d[p]
        d[p] = push!(v, l)
    else
        d[p] = [l]
    end

    return d
end


function PolACalculus(d::Dict{String,Array{Int16,1}})
    dframe = DataFrame(Chrmosome=String[], Position=Int32[], Median=Float16[], Minimum=Int16[], Maximum=Int16[], Counts=Int32[])

    for k in keys(d)
        md = median(d[k])
        mn = minimum(d[k])
        mx = maximum(d[k])
        l = length(d[k])
        push!(dframe, [split(k,"::")[1] parse(Int32,split(k,"::")[2]) md mn mx l])
    end

    return dframe
end


function WrFrame(fname::String, data::DataFrame)
    if !isfile(fname)
        CSV.write(fname, data)
    else
        CSV.write(fname, data, append=true, header=false)
    end
end


"Parse GFF file end extracts all exonic intervals from a gene - key of the returned dictionary - gene id"
function ParseGFF3(gff3file::String)::Dict{String,IntervalCollection{Bio.Intervals.GFF3.Record}}
    genes = Dict{String,IntervalCollection{Bio.Intervals.GFF3.Record}}()
    id = ""
    geneid = ""
    astr = GetAttributes(gff3file)

    for record in open(GFF3.Reader, gff3file)
        id = GFF3.attributes(record,"ID")[1]
        geneid = GFF3.attributes(record,"gene_id")[1]

        if  (GFF3.featuretype(record) == "exon" && String[id] in GFF3.attributes(record, "Parent"))
            if !(geneid in keys(genes))
                genes[geneid] = IntervalCollection{Bio.Intervals.GFF3.Record}()
            end

            push!(genes[geneid], Interval(record))
        end
    end

    return genes
end


function GetAttributes(gff3file::String)

    ct = Int16(0)
    s = Array{String,1}()

    for record in open(GFF3.Reader, gff3file)
        while ct < 10000

            for i in GFF3.attributes(record)
                if first(i) in s
                    continue
                else
                    push!(s, first(i))
                end
            end
            ct += 1
        end
    end

    return s
end


function main(args)
    arg_parse_settings = ArgParseSettings(description="Calculates uncovered windows number and GC bias.")
    @add_arg_table arg_parse_settings begin
        "--bam", "-b"
            nargs = '+'
            arg_type = String
            help = "Bam files separated by spaces."
            required = true
            dest_name = "bams"
        "--out", "-o"
            arg_type = String
            help = "Prefix of output file."
            required = false
            dest_name = "outfile"
            default = "output"
        "--gff", "-g"
            arg_type = String
            help = "GFF3 file"
            required = true
            dest_name = "gff3file"
    end

    parsed_args = parse_args(arg_parse_settings)

    for bam in parsed_args["bams"]
        treads, passreads, pareads, pclus = BamRead(bam)
        statdframe = PolACalculus(pclus)
        statdframe[:Smaple] = bam
        statdframe[:TotalReads] = treads
        statdframe[:PassedReads] = passreads
        statdframe[:PolyAReads] = pareads
        WrFrame(parsed_args["outfile"]*".csv", statdframe)
        genes = ParseGFF3(parsed_args["gff3file"])
    end
end

main(ARGS)
